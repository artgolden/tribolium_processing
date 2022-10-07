#@ String views_paths
#@ String psfs_paths
#@ String fused_output_path
#@ String num_deconv_iter
#@ String big_image_needs_batching
#@ ConvertService convert
#@ OpService ops
#@ DatasetService ds
#@ ConvertService convert

#############################################################################################################################################################
#       BEWARE! THIS SCRIPT LEAKS MEMORY! Like a lot... Intended to be run as a separate process, in its own headless Fiji instance.                        #
#############################################################################################################################################################

views_paths = views_paths
psfs_paths = psfs_paths
fused_output_path = fused_output_path
num_deconv_iter = int(num_deconv_iter)
ops = ops
unicode = unicode
ds = ds
convert = convert
big_image_needs_batching = big_image_needs_batching

# TODO: 
# - Check where 32-bit conversion is not needed
# - find CLIJ conversion to float


from datetime import datetime
import timeit
import os
import logging
from ij import IJ, ImagePlus
from ij.io import FileSaver
from net.haesleinhuepf.clijx import CLIJx
from net.haesleinhuepf.clij.clearcl import ClearCLBuffer
from net.haesleinhuepf.clijx.plugins import DeconvolveRichardsonLucyFFT
from net.haesleinhuepf.clij.coremem.enums import NativeTypeEnum


def logging_broadcast(string):
	print(string)
	logging.info(string)


def save_tiff_simple(image, path):
    if os.path.exists(path):
        os.remove(path)
    fs = FileSaver(image)
    fs.saveAsTiff(path)

def deconvolve_fuse_timepoint_multiview_entropy_weighted(views, transformed_psfs, num_iterations=15, sigma_scaling_factor_xy=1, sigma_scaling_factor_z=0.5, big_image_needs_batching=False):
    """Do deconvolution on individual views, then fuse them adjusting for entropy

    Args:
        views (ImagePlus[]): list of registered(transformed to be aligned) image stacks, ready for fusion
        transformed_psfs (ImagePlus[]): 32-bit transformed PSF for each view according to each view's registration affine transformation
        num_iterations (int, optional): number of deconvolution iterations. Defaults to 8.
        sigma_scaling_factor_xy (int, optional): scaling factor of X and Y axis, used to calculate sigmas for quick entropy calculation. Defaults to 4.
        sigma_scaling_factor_z (int, optional): scaling factor of Z axis, used to calculate sigmas for quick entropy calculation. Defaults to 2.
        big_image_needs_batching (bool, optional): whether to batch the deconvolution of the image. Defaults to False.

    Returns:
        ImagePlus: 32-bit fused image 
    """
    n_views = len(views)
    for i, view in enumerate(views):
        view_converted = ops.run("convert.float32", view)
        views[i] = view_converted

    start_time = timeit.default_timer()

    # init GPU
    clijx = CLIJx.getInstance()

    avg_weighted_image = clijx.convert(views[0], ClearCLBuffer)
    clijx.set(avg_weighted_image, 0)
    transformed_psfs_gpu = [clijx.convert(psf, ClearCLBuffer) for psf in transformed_psfs]
    buffer1 = clijx.create(avg_weighted_image)
    buffer2 = clijx.create(avg_weighted_image)
    avg_weights = clijx.create(avg_weighted_image)


    def compute_entropy_weight(view, output, sigma1=20, sigma2=40):
        clijx.gaussianBlur3D(view, buffer1, sigma1 * sigma_scaling_factor_xy, sigma1 * sigma_scaling_factor_xy, sigma1 * sigma_scaling_factor_z)
        clijx.addImagesWeighted(view, buffer1, buffer2, 1, -1)
        # clijx.convertFloat(buffer2, buffer1) # Is this needed? (the method does not exist :( ))
        clijx.multiplyImages(buffer2, buffer2, buffer1)
        clijx.gaussianBlur3D(buffer1, output, sigma2 * sigma_scaling_factor_xy, sigma2 * sigma_scaling_factor_xy, sigma2 * sigma_scaling_factor_z)

    # Weight views
    for i, psf in zip(range(n_views), transformed_psfs_gpu):
        logging_broadcast("Deconvolving %s view" % i)
        view = clijx.convert(views[i], ClearCLBuffer)
        clijx.release(buffer2)
        if(big_image_needs_batching):
            logging_broadcast("Recieved big_image_needs_batching parameter. Splitting in two batches for deconvolution")
            clijx.release(buffer2)
            avg_weighted_RAM = clijx.pull(avg_weighted_image)
            clijx.release(avg_weighted_image)
            avg_weights_RAM = clijx.pull(avg_weights)
            clijx.release(avg_weights)

            width, height, depth  = clijx.getDimensions(view)
            psf_width, psf_height, psf_depth = clijx.getDimensions(psf)
            half_height = height / 2 + psf_height / 2 + 2

            half_buffer_in = clijx.create([width, half_height, depth], NativeTypeEnum.Float)
            half_buffer_out = clijx.create(half_buffer_in)

            clijx.paste3D(view, half_buffer_in, 0, 0, 0)
            
            logging_broadcast("Deconvolving first half")
            DeconvolveRichardsonLucyFFT.deconvolveRichardsonLucyFFT(clijx, half_buffer_in, psf, half_buffer_out, num_iterations, 0.0, False)

            clijx.paste3D(view, half_buffer_in, 0, -(height - half_height), 0)
            clijx.paste3D(half_buffer_out, view, 0, 0, 0)
            
            logging_broadcast("Deconvolving second half")
            DeconvolveRichardsonLucyFFT.deconvolveRichardsonLucyFFT(clijx, half_buffer_in, psf, half_buffer_out, num_iterations, 0.0, False)

            clijx.paste3D(half_buffer_out, view, 0, height - half_height, 0)

            avg_weighted_image = clijx.convert(avg_weighted_RAM, ClearCLBuffer)
            avg_weights = clijx.convert(avg_weights_RAM, ClearCLBuffer)
        else:
            DeconvolveRichardsonLucyFFT.deconvolveRichardsonLucyFFT(clijx, view, psf, buffer1, num_iterations, 0.0, False)
            clijx.copy(buffer1, view)
        buffer2 = clijx.create(view)
        weight = buffer2
        compute_entropy_weight(view, weight)
        clijx.multiplyImages(view, weight, buffer1)

        # divide view for averaging
        clijx.addImagesWeighted(buffer1, avg_weighted_image, view, 1.0 / n_views, 1)
        clijx.copy(view, avg_weighted_image)

        clijx.addImagesWeighted(weight, avg_weights, buffer1, 1.0 / n_views, 1)
        clijx.copy(buffer1, avg_weights)
        clijx.release(view)


    # Normalize by average weights
    output = buffer2
    clijx.divideImages(avg_weighted_image, avg_weights, output)
    

    deconv_fused_image = clijx.pull(output)

    logging_broadcast("Deconvolution+content-based fusion took: %s ms" % round((timeit.default_timer() - start_time) * 1000, 1))
    logging_broadcast(clijx.reportMemory())
    # clean up
    clijx.clear()

    return deconv_fused_image


views_dir = os.path.dirname(views_paths.split(";")[0])



def get_dt_string():
    now = datetime.now()
    return now.strftime("%Y-%b-%d-%H%M%S")

dt_string = get_dt_string()
logging.basicConfig(filename=os.path.join(views_dir, "%s-TribFuse_CLIJ_deconv_fusion.log" % dt_string),
                    filemode='w',
                    format='%(asctime)s-%(levelname)s - %(message)s',
                    datefmt='%d-%b-%y %H:%M:%S',
                    level=logging.INFO)

logging_broadcast("Loading views")
start_time = timeit.default_timer()
views = []
for view_path in views_paths.split(";"):
    views.append(IJ.openImage(unicode(view_path)))
logging_broadcast("Loading views took: %s s" % round((timeit.default_timer() - start_time), 3))

logging_broadcast("Loading psfs")
transformed_psfs = []
for psf_path in psfs_paths.split(";"):
    transformed_psfs.append(IJ.openImage(unicode(psf_path)))

big_image_needs_batching = big_image_needs_batching == "True"

logging_broadcast("Starting deconvolution+fusion")
deconv_fused_image = deconvolve_fuse_timepoint_multiview_entropy_weighted(views, transformed_psfs, num_iterations=num_deconv_iter, big_image_needs_batching=big_image_needs_batching)

logging_broadcast("Converting output to 16-bit")
deconv_fused_image = ds.create(ops.run("convert.uint16", deconv_fused_image))
fused_imp = convert.convert(deconv_fused_image, ImagePlus)

logging_broadcast("Saving the output image")
save_tiff_simple(fused_imp, unicode(fused_output_path))

