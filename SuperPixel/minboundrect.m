function [rectx,recty,area,perimeter] = minboundrect(x,y,metric)
% minboundrect: Compute the minimal bounding rectangle of points in the plane计算平面中点的最小边界矩形
% usage: [rectx,recty,area,perimeter] = minboundrect(x,y,metric)
%
% arguments: (input)
%  x,y - vectors of points, describing points in the plane as      点的矢量，将平面中的点描述
%        (x,y) pairs. x and y must be the same lengths.为（x，y）   对。 x和y的长度必须相同。
%
%  metric - (OPTIONAL) - single letter character flag which   （可选） - 单字母字符标志，  
%        denotes the use of minimal area or perimeter as the   表示使用最小面积或周长作为最小化
%        metric to be minimized. metric may be either 'a' or 'p',  的度量。度量可以是'a'或'p'，
%        capitalization is ignored. Any other contraction of 'area'  大写被忽略。
%        or 'perimeter' is also accepted.   “区域”或“周长”的任何其他收缩也被接受。
%
%        DEFAULT: 'a'    ('area')
%
% arguments: (output)
%  rectx,recty - 5x1 vectors of points that define the minimal
%        bounding rectangle.  定义最小边界矩形的5x1点矢量。
%
%  area - (scalar) area of the minimal rect itself. （标量）最小矩形本身的区域。
%
%  perimeter - (scalar) perimeter of the minimal rect as found 找到的最小矩形的（标量）周长
%
%
% Note: For those individuals who would prefer the rect with minimum   对于那些喜欢具有最小周长或  
% perimeter or area, careful testing convinces me that the minimum area  面积的矩形的人来说，仔细的
% rect was generally also the minimum perimeter rect on most problems  测试使我确信最小区域rect通常也
% (with one class of exceptions). This same testing appeared to verify my   是大多数问题的最小周长（有
% assumption that the minimum area rect must always contain at least    一类例外）。同样的测试似乎证实了
% one edge of the convex hull. The exception I refer to above is for    我的假设，即最小区域rect必须始终
% problems when the convex hull is composed of only a few points,    包含凸包的至少一个边缘。我在上面提到
% most likely exactly 3. Here one may see differences between the    的例外是当凸包仅由几个点组成时的问题，
% two metrics. My thanks to Roger Stafford for pointing out this     很可能恰好是3.这里可以看到两个指标
% class of counter-examples.   之间的差异。 我要感谢Roger Stafford指出这类反例。
%
% Thanks are also due to Roger for pointing out a proof that the谢谢也是由于Roger在
% bounding rect must always contain an edge of the convex hull, in最小周长和区域情况下指出了
% both the minimal perimeter and area cases.边界矩形必须始终包含凸包边缘的证明。
%
%
% See also: minboundcircle, minboundtri, minboundsphere
%
%
% default for metric指标的默认值
if (nargin<3) || isempty(metric)  %nargin参数
  metric = 'a';
elseif ~ischar(metric)
  error 'metric must be a character flag if it is supplied.'
else
  % check for 'a' or 'p'
  metric = lower(metric(:)');  %小写                  
  ind = strmatch(metric,{'area','perimeter'});             
  if isempty(ind)                
    error 'metric does not match either ''area'' or ''perimeter'''
  end
  % just keep the first letter.
  metric = metric(1);
end
 
 
% preprocess data
x=x(:);
y=y(:);
 
 
% not many error checks to worry about
n = length(x);                                    
if n~=length(y)                               
  error 'x and y must be the same sizes'
end
 
 
% start out with the convex hull of the points to  从点的凸包开始，以显着减少
% reduce the problem dramatically. Note that any  问题。 请注意，凸包内部的
% points in the interior of the convex hull are   任何点都不需要，
% never needed, so we drop them.    所以我们放下它们。
if n>3
  edges = convhull(x,y);  % 'Pp' will silence the warnings
                          %计算轮廓凸包
 
 
  % exclude those points inside the hull as not relevant
  % also sorts the points into their convex hull as a
  % closed polygon
 
  x = x(edges);
  y = y(edges);
 
  % probably fewer points now, unless the points are fully convex
  nedges = length(x) - 1;                       
elseif n>1
  % n must be 2 or 3
  nedges = n;
  x(end+1) = x(1);
  y(end+1) = y(1);
else
  % n must be 0 or 1
  nedges = n;
end
 
 
% now we must find the bounding rectangle of those
% that remain.
 
 
% special case small numbers of points. If we trip any
% of these cases, then we are done, so return.
switch nedges
  case 0
    % empty begets empty
    rectx = [];
    recty = [];
    area = [];
    perimeter = [];
    return
  case 1
    % with one point, the rect is simple.
    rectx = repmat(x,1,5);
    recty = repmat(y,1,5);
    area = 0;
    perimeter = 0;
    return
  case 2
    % only two points. also simple.
    rectx = x([1 2 2 1 1]);
    recty = y([1 2 2 1 1]);
    area = 0;
    perimeter = 2*sqrt(diff(x).^2 + diff(y).^2);
    return
end
% 3 or more points.
 
 
% will need a 2x2 rotation matrix through an angle theta
Rmat = @(theta) [cos(theta) sin(theta);-sin(theta) cos(theta)];
 
 
% get the angle of each edge of the hull polygon.
ind = 1:(length(x)-1);
edgeangles = atan2(y(ind+1) - y(ind),x(ind+1) - x(ind));
% move the angle into the first quadrant.
edgeangles = unique(mod(edgeangles,pi/2));
 
 
% now just check each edge of the hull
nang = length(edgeangles);              
area = inf;                           
perimeter = inf;
met = inf;
xy = [x,y];
for i = 1:nang                         
  % rotate the data through -theta 
  rot = Rmat(-edgeangles(i));
  xyr = xy*rot;
  xymin = min(xyr,[],1);
  xymax = max(xyr,[],1);
 
  % The area is simple, as is the perimeter
  A_i = prod(xymax - xymin);
  P_i = 2*sum(xymax-xymin);
 
  if metric=='a'
    M_i = A_i;
  else
    M_i = P_i;
  end
 
  % new metric value for the current interval. Is it better?
  if M_i<met
    % keep this one
    met = M_i;
    area = A_i;
    perimeter = P_i;
 
    rect = [xymin;[xymax(1),xymin(2)];xymax;[xymin(1),xymax(2)];xymin];
    rect = rect*rot';
    rectx = rect(:,1);
    recty = rect(:,2);
  end
end