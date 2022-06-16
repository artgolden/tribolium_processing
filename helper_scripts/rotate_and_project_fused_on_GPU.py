#@ Dataset data
#@ Integer rotation_angleY
#@OUTPUT projections_stack
# OpService ops
# DatasetService ds
# ConvertService convert
# UIService ui


import timeit
from ij import IJ
from net.haesleinhuepf.clij2 import CLIJ2
from net.haesleinhuepf.clij.clearcl import ClearCLBuffer



data = data
stack_dims = data.dimensionsAsLongArray()

start_time = timeit.default_timer()

# init GPU
clij2 = CLIJ2.getInstance()

# push image to GPU
input = clij2.convert(data, ClearCLBuffer)

rotated = clij2.create(input)
rotated_half = clij2.create([stack_dims[0], stack_dims[1], stack_dims[2] / 2], input.getNativeType())
destination_max = clij2.create(rotated_half)
projection_2d = clij2.create([stack_dims[0], stack_dims[1]], input.getNativeType())

angleX = 0.0
angleY = rotation_angleY
angleZ = 0
rotateAroundCenter = True
num_projections = 360/int(angleY)

projections = clij2.create([stack_dims[0], stack_dims[1], num_projections], input.getNativeType())

for i in range(num_projections):
	angleY = rotation_angleY * (i + 1)
	clij2.rotate3D(input, rotated, angleX, angleY, angleZ, rotateAroundCenter)
	clij2.copy( rotated, rotated_half)

	clij2.maximumZProjection(rotated_half, destination_max);
	clij2.copy(destination_max, projection_2d)
	clij2.paste(projection_2d, projections, 0, 0, i)


projections_stack = clij2.pull(projections)


print("Total execution time: %s ms" % round((timeit.default_timer() - start_time) * 1000, 1))
print(clij2.reportMemory())
# clean up
clij2.clear()