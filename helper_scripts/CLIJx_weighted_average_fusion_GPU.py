# Dataset data
# Dataset psf
# Integer num_iterations
#OUTPUT projections_stack


# TODO: 
# - Check where 32-bit conversion is not needed
# - find CLIJ conversion to float

import timeit
from ij import IJ
from net.haesleinhuepf.clijx import CLIJx
from net.haesleinhuepf.clij.clearcl import ClearCLBuffer
from net.haesleinhuepf.clijx.plugins import DeconvolveRichardsonLucyFFT

IJ.run("Close All")
direction1_raw = IJ.openImage("/media/tema/big_storage/work/goethe/fiji/test_dataset/DS0016_MEME6/downsampled_23_timepoint/raw_cropped_registered_for_fusion/fused_tp_0_vs_0.tif")
direction2_raw = IJ.openImage("/media/tema/big_storage/work/goethe/fiji/test_dataset/DS0016_MEME6/downsampled_23_timepoint/raw_cropped_registered_for_fusion/fused_tp_0_vs_1.tif")
direction3_raw = IJ.openImage("/media/tema/big_storage/work/goethe/fiji/test_dataset/DS0016_MEME6/downsampled_23_timepoint/raw_cropped_registered_for_fusion/fused_tp_0_vs_2.tif")
direction4_raw = IJ.openImage("/media/tema/big_storage/work/goethe/fiji/test_dataset/DS0016_MEME6/downsampled_23_timepoint/raw_cropped_registered_for_fusion/fused_tp_0_vs_3.tif")

dataset_xml_path = "/media/tema/big_storage/work/goethe/fiji/test_dataset/DS0016_MEME6/downsampled_23_timepoint/dataset.xml"

transformed_psf_view_1 = IJ.openImage("/media/tema/big_storage/work/goethe/fiji/test_dataset/DS0016_MEME6/downsampled_23_timepoint/transformed_psf/tp_0_view_1.tif")
transformed_psf_view_2 = IJ.openImage("/media/tema/big_storage/work/goethe/fiji/test_dataset/DS0016_MEME6/downsampled_23_timepoint/transformed_psf/tp_0_view_2.tif")
transformed_psf_view_3 = IJ.openImage("/media/tema/big_storage/work/goethe/fiji/test_dataset/DS0016_MEME6/downsampled_23_timepoint/transformed_psf/tp_0_view_3.tif")
transformed_psf_view_4 = IJ.openImage("/media/tema/big_storage/work/goethe/fiji/test_dataset/DS0016_MEME6/downsampled_23_timepoint/transformed_psf/tp_0_view_4.tif")

reference_view = direction1_raw

views = [direction1_raw, direction2_raw, direction3_raw, direction4_raw]
transformed_psfs = [transformed_psf_view_1, transformed_psf_view_2, transformed_psf_view_3, transformed_psf_view_4]

for imp in views:
    IJ.run(imp, "32-bit", "")

n_views = len(views)

downsampling_factor_xy = 4
num_iterations = 8

start_time = timeit.default_timer()

# init GPU
clijx = CLIJx.getInstance()

# views_gpu = [clijx.convert(view, ClearCLBuffer) for view in views]
avg_weighted_image = clijx.convert(views[0], ClearCLBuffer)
clijx.set(avg_weighted_image, 0)
transformed_psfs_gpu = [clijx.convert(psf, ClearCLBuffer) for psf in transformed_psfs]
buffer1 = clijx.create(avg_weighted_image)
buffer2 = clijx.create(avg_weighted_image)
avg_weights = clijx.create(avg_weighted_image)




def compute_entropy_weight(view, output, sigma1=20, sigma2=40):
    clijx.gaussianBlur3D(view, buffer1, sigma1 / downsampling_factor_xy, sigma1 / downsampling_factor_xy, sigma1 / 2)
    clijx.addImagesWeighted(view, buffer1, buffer2, 1, -1)
    # clijx.convertFloat(buffer2, buffer1) # Is this needed? (the method does not exist :( ))
    clijx.multiplyImages(buffer2, buffer2, buffer1)
    clijx.gaussianBlur3D(buffer1, output, sigma2 / downsampling_factor_xy, sigma2 / downsampling_factor_xy, sigma2 / 2)

# Weight views
for i, psf in zip(range(n_views), transformed_psfs_gpu):
    print("Deconvolving %s view" % i)
    view = clijx.convert(views[i], ClearCLBuffer)
    DeconvolveRichardsonLucyFFT.deconvolveRichardsonLucyFFT(clijx, view, psf, buffer1, num_iterations, 0.0, False)
    clijx.copy(buffer1, view)
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



    

clijx.pull(output).show()

print("Total execution time: %s ms" % round((timeit.default_timer() - start_time) * 1000, 1))
print(clijx.reportMemory())
# clean up
clijx.clear()