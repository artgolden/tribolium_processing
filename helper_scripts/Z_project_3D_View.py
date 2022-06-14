#@ Dataset data
#@ String(label="Projection Type", choices={"Max","Mean","Median","Min", "StdDev", "Sum"}) projection_type
#@OUTPUT Dataset output
#@ OpService ops
#@ DatasetService ds

# Project 3D View stack along Z dimension

from net.imagej.axis import Axes
from net.imagej.ops import Ops

# Select which dimension to project
dim = 2


# Write the output dimensions
new_dimensions = [data.dimension(d) for d in range(0, data.numDimensions()) if d != dim]

# Create the output image
projected = ops.create().img(new_dimensions)

# Create the op and run it
proj_op = ops.op(getattr(Ops.Stats, projection_type), data)
ops.transform().project(projected, data, proj_op, dim)

# Create the output Dataset
output = ds.create(projected)
