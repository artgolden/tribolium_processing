#@ File (label='Image to rotate and project ') image_path
#@ File(label='Directory to save images', style='directory') output_dir
#@ Double (value=22.5, stepSize=0.5, label='Angle of rotation around Y-axis in degrees', style="format:#.00") rotation_angleY 
#@ Double (value=0, stepSize=0.5, label='Additional pre-rotation around Y-axis in degrees', style="format:#.00") rotation_addition 
#@ OpService ops
#@ DatasetIOService io
#@ DatasetService ds
# ConvertService convert
#@ UIService ui


import os
import timeit
import math
from ij import IJ
from ij.io import FileSaver
from net.haesleinhuepf.clij2 import CLIJ2
from net.haesleinhuepf.clij.clearcl import ClearCLBuffer

def save_tiff_simple(image, path):
    if os.path.exists(path):
        os.remove(path)
    fs = FileSaver(image)
    fs.saveAsTiff(path)



ds = ds
io = io
ops = ops
ui = ui
rotation_angleY = rotation_angleY
rotation_addition = rotation_addition
print("Got rotation angle %s" % rotation_angleY)
output_dir = output_dir.getAbsolutePath()
image_path = image_path.getAbsolutePath()

print("Loading the image %s" % image_path)
data = io.open(image_path)

data32bit = ops.run("convert.float32", data)
stack_dims = data32bit.dimensionsAsLongArray()

start_time = timeit.default_timer()

# init GPU
clij2 = CLIJ2.getInstance()

print("Creating images on the GPU")
# push image to GPU
input = clij2.convert(data32bit, ClearCLBuffer)

rotated = clij2.create(input)
destination_max = clij2.create(rotated)
rotated_adjusted = clij2.create(rotated)
projection_2d = clij2.create([stack_dims[0], stack_dims[1]], input.getNativeType())
depth_adjustment_stack = clij2.create([stack_dims[0], stack_dims[1], stack_dims[2]])
adjustment_2d = projection_2d

print("Creating adjustment stack")
nplanes = stack_dims[2]
for plane in range(nplanes):
	clij2.set(adjustment_2d, float(nplanes - plane) / nplanes)

	clij2.paste(adjustment_2d, depth_adjustment_stack, 0, 0, plane)

angleX = 0
angleY = rotation_angleY
angleZ = 0
rotateAroundCenter = True
num_projections = 360/int(angleY)

print("Rotating and creating projections")
projections = clij2.create([stack_dims[0], stack_dims[1], num_projections], input.getNativeType())
for i in range(num_projections):
	angleY = (rotation_angleY * i + rotation_addition ) * math.pi / 180  # you need to switch to angles for the newer version of CLIJ
	print("Rotating by angle: %s" % (rotation_angleY * i))
	clij2.rotate3D(input, rotated, angleX, angleY, angleZ, rotateAroundCenter)
	clij2.multiplyImages(rotated, depth_adjustment_stack, rotated_adjusted)
	clij2.maximumZProjection(rotated_adjusted, destination_max);
	clij2.copy(destination_max, projection_2d)
	clij2.paste(projection_2d, projections, 0, 0, i)


projections_stack = clij2.pull(projections)
# print (projections_stack, projections_stack.getClass())
projections_stack = ds.create((ops.run("convert.uint16", projections_stack)))

print("Total execution time: %s ms" % round((timeit.default_timer() - start_time) * 1000, 1))
print(clij2.reportMemory())
# clean up
clij2.clear()
output_path = os.path.join(output_dir, os.path.basename(image_path))
if os.path.exists(output_path):
	os.remove(output_path)
	print("removing existing image")
print("Saving projections")
io.save(projections_stack, output_path)