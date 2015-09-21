import math;

def webmercator2wgs84(x):
	mercatorX_lon=x[0];
	mercatorY_lat=x[1];
	if (abs(mercatorX_lon) < 180 and abs(mercatorY_lat) < 90):
		return;
	if ((abs(mercatorX_lon) > 20037508.3427892) or (abs(mercatorY_lat) > 20037508.3427892)):
		return;
	x = mercatorX_lon;
	y = mercatorY_lat;
	num3 = x / 6378137.0;
	num4 = num3 * 57.295779513082323;
	num5 = math.floor(((num4 + 180.0) / 360.0));
	num6 = num4 - (num5 * 360.0);
	num7 = 1.5707963267948966 - (2.0 * math.atan(math.exp((-1.0 * y) / 6378137.0)));
	mercatorX_lon = num6;
	mercatorY_lat = num7 * 57.295779513082323;
	return [mercatorX_lon,mercatorY_lat];


#private void ToWebMercator(ref double mercatorX_lon, ref double mercatorY_lat)
#{
#    if ((abs(mercatorX_lon) > 180 || abs(mercatorY_lat) > 90))
#        return;#

    #double num = mercatorX_lon * 0.017453292519943295;
    #double x = 6378137.0 * num;
    #double a = mercatorY_lat * 0.017453292519943295;

   # mercatorX_lon = x;#
   # mercatorY_lat = 3189068.5 * Math.Log((1.0 + Math.Sin(a)) / (1.0 - Math.Sin(a)));
#}

#near munich
#x=[1252387.6317574722, 6102855.354436191]
#print x;
#y = webmercator2wgs84(x);
#print y;
