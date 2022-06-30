#@ Dataset data
#@ Integer rotation_angleY
#@OUTPUT projections_stack
#@ OpService ops
# DatasetService ds
# ConvertService convert
#@ UIService ui


import timeit
import math
from ij import IJ
from net.haesleinhuepf.clij2 import CLIJ2
from net.haesleinhuepf.clij.clearcl import ClearCLBuffer



data = data
data32bit = ops.run("convert.float32", data)
stack_dims = data32bit.dimensionsAsLongArray()

start_time = timeit.default_timer()

# init GPU
clij2 = CLIJ2.getInstance()

# push image to GPU
input = clij2.convert(data32bit, ClearCLBuffer)

rotated = clij2.create(input)
destination_max = clij2.create(rotated)
rotated_adjusted = clij2.create(rotated)
projection_2d = clij2.create([stack_dims[0], stack_dims[1]], input.getNativeType())
depth_adjustment_stack = clij2.create([stack_dims[0], stack_dims[1], stack_dims[2]])
adjustment_2d = clij2.create([stack_dims[0], stack_dims[1]], input.getNativeType())

nplanes = stack_dims[2]
for plane in range(nplanes):
	clij2.set(adjustment_2d, float(nplanes - plane) / nplanes)

	clij2.paste(adjustment_2d, depth_adjustment_stack, 0, 0, plane)

angleX = 0.0
angleY = rotation_angleY
angleZ = 0
rotateAroundCenter = True
num_projections = 360/int(angleY)

projections = clij2.create([stack_dims[0], stack_dims[1], num_projections], input.getNativeType())

for i in range(num_projections):
	angleY = rotation_angleY * i * math.pi / 180
	clij2.rotate3D(input, rotated, angleX, angleY, angleZ, rotateAroundCenter)
	clij2.multiplyImages(rotated, depth_adjustment_stack, rotated_adjusted)
	clij2.maximumZProjection(rotated_adjusted, destination_max);
	clij2.copy(destination_max, projection_2d)
	clij2.paste(projection_2d, projections, 0, 0, i)


projections_stack = clij2.pull(projections)
projections_stack = ops.run("convert.uint16", projections_stack)


print("Total execution time: %s ms" % round((timeit.default_timer() - start_time) * 1000, 1))
print(clij2.reportMemory())
# clean up
clij2.clear()