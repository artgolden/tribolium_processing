#@ Dataset data
#@ OpService ops
#@ DatasetService ds
#@ ConvertService convert
#@ UIService ui

from numbers import Real
import jarray
from net.imglib2 import Cursor
from net.imglib2 import RandomAccessible
from net.imglib2 import RandomAccessibleInterval
from net.imglib2 import RealRandomAccessible
from net.imglib2.img import Img
from net.imglib2.img.array import ArrayImgFactory
from net.imglib2.img.display.imagej import ImageJFunctions
from net.imglib2.interpolation.randomaccess import NLinearInterpolatorFactory
from net.imglib2.realtransform import AffineTransform
from net.imglib2.realtransform import RealViews
from net.imglib2.type.numeric import NumericType
from net.imglib2.type.numeric.real import FloatType
from net.imglib2.view import Views
from net.imglib2 import FinalInterval
from net.imglib2.interpolation.randomaccess import NLinearInterpolatorFactory
from net.imglib2.realtransform import AffineTransform3D
from net.imglib2.realtransform import RealViews
from net.imglib2.util import Intervals
from net.imglib2.view import Views
from net.imagej.ops.transform.translateView import IntervalTranslateView

x_size = 521
y_size = 1025
z_size = 521

region = FinalInterval.createMinSize( 0, 0, 0,    x_size, y_size, z_size)

affine_to_center = AffineTransform3D()
affine_to_center.translate( -x_size/2, -y_size/2, -z_size/2 );

affine_rotate = AffineTransform3D()
affine_rotate.rotate(1, (3.14159 / 4) );

affine_back = AffineTransform3D()
affine_back.translate( x_size/2, y_size/2, z_size/2 );

interpolator = NLinearInterpolatorFactory()
# translated = ops.op("translateView", data, jarray.array([-60, -123, -525], 'l'))
# translated = IntervalTranslateView(data, jarray.array([-60, -123, -525], 'l'))
# ui.show(translated);
# exit()

ext_view = Views.extendZero( data ); # extend
# ext_view = ops.run("translateView", ext_view)
# ext_view = Views.translate(ext_view, [-60, -123, -525])
real_view = Views.interpolate( ext_view, interpolator );
affine_box = RealViews.affine( real_view, affine_to_center );
affine_box = RealViews.affine( affine_box, affine_rotate );
affine_box = RealViews.affine( affine_box, affine_back );
affine_view = Views.interval( affine_box, region );
ui.show(affine_view);
