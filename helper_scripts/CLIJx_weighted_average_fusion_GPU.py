# Dataset data
# Dataset psf
# Integer num_iterations
#OUTPUT projections_stack


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

reference_view = direction1_raw

views = [direction1_raw, direction2_raw, direction3_raw, direction4_raw]
for imp in views:
    IJ.run(imp, "32-bit", "")

n_views = len(views)

downsampling_factor = 4


start_time = timeit.default_timer()

# init GPU
clijx = CLIJx.getInstance()

views_gpu = [clijx.convert(view, ClearCLBuffer) for view in views]
buffer1 = clijx.create(views_gpu[0])
buffer2 = clijx.create(views_gpu[0])
avg_weights = clijx.create(views_gpu[0])

# input_psf = clijx.convert(psf, ClearCLBuffer)
# DeconvolveRichardsonLucyFFT.deconvolveRichardsonLucyFFT(clijx, input, input_psf, output_deconv, num_iterations, 0.0, False)


def compute_entropy_weight(view, output, sigma1=20, sigma2=40):
    clijx.gaussianBlur3D(view, buffer1, sigma1 / downsampling_factor, sigma1 / downsampling_factor, sigma1 / 2)
    clijx.addImagesWeighted(view, buffer1, buffer2, 1, -1)
    # clijx.convertFloat(buffer2, buffer1) # Is this needed? (the method does not exist :( ))
    clijx.multiplyImages(buffer2, buffer2, buffer1)
    clijx.gaussianBlur3D(buffer1, output, sigma2 / downsampling_factor, sigma2 / downsampling_factor, sigma2 / 2)

# Weight views
for view in views_gpu:
    weight = buffer2
    compute_entropy_weight(view, weight)
    clijx.multiplyImages(view, weight, buffer1)

    # divide view for averaging
    clijx.multiplyImageAndScalar(buffer1, view, 1.0 / n_views)

    clijx.addImagesWeighted(weight, avg_weights, buffer1, 1.0 / n_views, 1)
    clijx.copy(buffer1, avg_weights)

# add divided weighted images to average
avg_weighted_img = buffer1
clijx.addImages(views_gpu[0], views_gpu[1], avg_weighted_img)
for i in range(2, n_views):
    clijx.addImages(views_gpu[i], avg_weighted_img, buffer2)
    clijx.copy(buffer2, avg_weighted_img)

# Normalize by average weights
output = buffer2
buffer4 = views_gpu[0]
clijx.divideImages(avg_weighted_img, avg_weights, output)



    

clijx.pull(output).show()

print("Total execution time: %s ms" % round((timeit.default_timer() - start_time) * 1000, 1))
print(clijx.reportMemory())
# clean up
clijx.clear()