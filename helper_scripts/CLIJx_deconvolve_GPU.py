#@ File(label='Choose an image to deconvolve', style='file') data_file
#@ File(label='Choose a PSF to deconvolve with', style='file') psf_file
#@ File(label='Choose a directory to save deconvolved images', style='directory') save_dir
#@ Integer num_iterations
#@ DatasetIOService io
#@ ConvertService convert


import timeit
import os
from ij import IJ
from net.haesleinhuepf.clijx import CLIJx
from net.haesleinhuepf.clij.clearcl import ClearCLBuffer
from net.haesleinhuepf.clijx.plugins import DeconvolveRichardsonLucyFFT
from net.imagej import Dataset
from ij.io import FileSaver

def save_tiff(image, path):
	if os.path.exists(path):
		os.remove(path)
	fs = FileSaver(image)
	fs.saveAsTiff(path)


data = io.open(data_file.getAbsolutePath())
psf = io.open(psf_file.getAbsolutePath())
stack_dims = data.dimensionsAsLongArray()

start_time = timeit.default_timer()

# init GPU
clijx = CLIJx.getInstance()

# push image to GPU
input = clijx.convert(data, ClearCLBuffer)
input_psf = clijx.convert(psf, ClearCLBuffer)
output_deconv = clijx.create(input)

DeconvolveRichardsonLucyFFT.deconvolveRichardsonLucyFFT(clijx, input, input_psf, output_deconv, num_iterations, 0.0, False)

save_path = os.path.join(save_dir.getAbsolutePath(), os.path.basename(data_file.getAbsolutePath()))

deconv_image = clijx.pull(output_deconv)
# io.save(convert.convert(clijx.pull(output_deconv), Dataset), save_path)
save_tiff(deconv_image, save_path)

print("Total execution time: %s ms" % round((timeit.default_timer() - start_time) * 1000, 1))
print(clijx.reportMemory())
# clean up
clijx.clear()