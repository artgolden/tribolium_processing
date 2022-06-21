import timeit
from ij import IJ
from net.haesleinhuepf.clijx import CLIJx
from net.haesleinhuepf.clij.clearcl import ClearCLBuffer
from net.haesleinhuepf.clijx.plugins import DeconvolveRichardsonLucyFFT


def deconvolve_fuse_timepoint_multiview_entropy_weighted(views, transformed_psfs, num_iterations=8, sigma_scaling_factor_xy=1, sigma_scaling_factor_z=0.5):
    """Do deconvolution on individual views, then fuse them adjusting for entropy

    Args:
        views (ImagePlus[]): 32-bit list of registered(transformed to be aligned) image stacks, ready for fusion
        transformed_psfs (ImagePlus[]): 32-bit transformed PSF for each view according to each view's registration affine transformation
        num_iterations (int, optional): number of deconvolution iterations. Defaults to 8.
        sigma_scaling_factor_xy (int, optional): scaling factor of X and Y axis, used to calculate sigmas for quick entropy calculation. Defaults to 4.
        sigma_scaling_factor_z (int, optional): scaling factor of Z axis, used to calculate sigmas for quick entropy calculation. Defaults to 2.

    Returns:
        ImagePlus: 32-bit fused image 
    """
    n_views = len(views)

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
        print("Deconvolving %s view" % i)
        view = clijx.convert(views[i], ClearCLBuffer)
        clijx.release(buffer2)
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

    print("Deconvolution+content-based fusion took: %s ms" % round((timeit.default_timer() - start_time) * 1000, 1))
    print(clijx.reportMemory())
    # clean up
    clijx.clear()

    return deconv_fused_image