// ImagePlus imp

run("CLIJ2 Macro Extensions", "cl_device=[Quadro P1000]");

// maximum y projection
image1 = getTitle();
Ext.CLIJ2_push(image1);
image2 = "maximum_y_projection-" + image1;
Ext.CLIJ2_maximumYProjection(image1, image2);
Ext.CLIJ2_pull(image2);
run("Canvas Size...", "width=640 height=200 position=Center zero");

// scale
Ext.CLIJ2_push(image2);
image3 = "scale-" + image1;
scaling_factor_x = 1.0;
scaling_factor_y = 4.0;
scale_to_center = true;
Ext.CLIJ2_scale2D(image2, image3, scaling_factor_x, scaling_factor_y, scale_to_center);
Ext.CLIJ2_pull(image3);
Ext.CLIJ2_clear();