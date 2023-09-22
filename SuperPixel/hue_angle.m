function h=hue_angle(a,b)
% HUE_ANGLE: Computes four-quadrant polar angle in degrees from Cartesian coordinates 
%
%   Colour Engineering Toolbox
%   author:    © Phil Green
%   version:   1.1
%   date:  	   17-01-2001
%   book:      http://www.wileyeurope.com/WileyCDA/WileyTitle/productCd-0471486884.html
%   web:       http://www.digitalcolour.org

h=(180/pi)*atan2(b,a);
j=(b<0);
h(j)=h(j)+360;