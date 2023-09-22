% CIRCULARSTRUCT
%
% Function to construct a circular structuring element
% for morphological operations.  构造用于形态学运算的圆形结构元素的功能。
%
% function strel = circularstruct(radius)
%
% Note radius can be a floating point value though the resulting
% circle will be a discrete approximation 注意半径可以是浮点值，但结果圆将是离散近似值
%
% Peter Kovesi   March 2000

function strel = circularstruct(radius)

if radius < 1
  error('radius must be >= 1');
end

dia = ceil(2*radius);  % Diameter of structuring element结构元素的直径
                        %ceil朝正无穷大方向取整

if mod(dia,2) == 0     % If diameter is a odd value如果直径是奇数值
 dia = dia + 1;        % add 1 to generate a `centre pixel'加1以生成“中心像素”
end

r = fix(dia/2);%向零方向取整
[x,y] = meshgrid(-r:r);
rad = sqrt(x.^2 + y.^2);  %x.^2每个元素平方
strel = rad <= radius;

