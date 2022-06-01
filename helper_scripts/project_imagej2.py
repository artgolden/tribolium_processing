#@ ImagePlus imp
#@ String(label="Dimension to Project", choices={"X", "Y", "Z", "TIME", "CHANNEL"}) projected_dimension
#@ String(label="Projection Type", choices={"Max","Mean","Median","Min", "StdDev", "Sum"}) projection_type
#@OUTPUT ImagePlus output
#@ OpService ops
#@ DatasetService ds
#@ ConvertService convert

# Do a projection of a stack. The projection is done along a specified axis.
# The commin use case of this script is to do a maximum Z projection.
 
from net.imagej.axis import Axes
from net.imagej.ops import Ops
from net.imagej import Dataset
from ij import IJ, ImagePlus, ImageStack, WindowManager

 

def project_image(image, dim_to_project, projection_type):
    data = convert.convert(image, Dataset)


    # Select which dimension to project
    dim = data.dimensionIndex(getattr(Axes, dim_to_project))

    if dim == -1:
        raise Exception("%s dimension not found." % dim_to_project)

    if data.dimension(dim) < 2:
        raise Exception("%s dimension has only one frame." % dim_to_project)

    # Write the output dimensions
    new_dimensions = [data.dimension(d) for d in range(0, data.numDimensions()) if d != dim]

    # Create the output image
    projected = ops.create().img(new_dimensions)


    # Create the op and run it
    proj_op = ops.op(getattr(Ops.Stats, projection_type), data)
    ops.transform().project(projected, data, proj_op, dim)

    projected = ops.run("convert.uint16", projected)

    out = ds.create(projected)
    o = convert.convert(out, ImagePlus)
    ip = o.getProcessor()
    o = ImagePlus("projection", ip)
    return o

output = project_image(imp, projected_dimension, projection_type)
