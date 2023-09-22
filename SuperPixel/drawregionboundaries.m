% DRAWREGIONBOUNDARIES Draw boundaries of labeled regions in an image     DRAWREGIONBOUNDARIES绘制图像中标记区域的边界
%
% Usage: maskim = drawregionboundaries(l, im, col)
%
% Arguments:参数
%            l - Labeled image of regions.区域标记的图像
%           im - Optional image to overlay the region boundaries on.可选图像以覆盖区域边界
%          col - Optional colour specification. Defaults to black.  Note that可选颜色规格。 默认为黑色。 请注意，
%                the colour components are specified as values 0-255.颜色分量指定为值0-255。
%                For example red is [255 0 0] and white is [255 255 255].例如，红色是[255 0 0]，白色是[255 255 255]。
%
% Returns: 
%       maskim - If no image has been supplied maskim is a binary mask如果没有提供图像，则maskim是二进制掩模图像，
%                image indicating where the region boundaries are.指示区域边界的位置。
%                If an image has been supplied maskim is the image with the如果提供了图像，则maskim是覆盖区域边界的图像
%                region boundaries overlaid 
%
% See also: MASKIMAGE

% Copyright (c) 2013 Peter Kovesi
% www.peterkovesi.com/matlabfns/
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in 
% all copies or substantial portions of the Software.
%
% The Software is provided "as is", without warranty of any kind.

% Feb 2013

function maskim = drawregionboundaries(l, im, col)
    
    % Form the mask by applying a sobel edge detector to the labeled image,通过将sobel边缘检测器应用于标记图像来形成掩模，
    % thresholding and then thinning the result.阈值处理然后稀疏结果。
%    h = [1  0 -1
%         2  0 -2
%         1  0 -1];
    h = [-1 1];  % A simple small filter is better in this application.
                 % Small regions 1 pixel wide get missed using a Sobel
                 % operator 一个简单的小滤波器在这个应用中更好。使用Sobel算子可以丢失1个像素宽的小区域
% h = fspecial('prewitt');%滤波算子
    gx = filter2(h ,l);
    gy = filter2(h',l);
    maskim = (gx.^2 + gy.^2) > 0;
    maskim = bwmorph(maskim, 'thin', Inf);%对二值图像maskim进行细化操作
    
    % Zero out any mask values that may have been set around the edge of the
    % image.  清零可能已在图像边缘周围设置的任何掩模值。
    maskim(1,:) = 0; maskim(end,:) = 0;
    maskim(:,1) = 0; maskim(:,end) = 0;
    
    % If an image has been supplied apply the mask to the image and return it 如果已提供图像，请将掩模应用于图像并将其返回
    if exist('im', 'var') 
        if ~exist('col', 'var'), col = 0; end
        maskim = maskimage(im, maskim, col);
    end