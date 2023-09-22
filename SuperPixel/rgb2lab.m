% RGB2LAB - RGB to L*a*b* colour space
%
% Usage: Lab = rgb2lab(im, wp)
%
% Arguments:  im - RGB image or Nx3 colourmap for conversion用于转换的RGB图像或Nx3色彩图
%             wp - Optional string specifying the adapted white point.可选字符串，指定适应的白点。
%                  This defaults to 'D65'.默认为“D65”
%
% Returns:   Lab - The converted image or colourmap.转换后的图像或颜色图
%
% This function wraps up calls to MAKECFORM and APPLYCFORM in a convenient此函数以方便的形式包含对MAKECFORM和APPLYCFORM的调用。 
% form.  Note that if the image is of type uint8 this function casts it to请注意，如果图像的类型为uint8，
% double and divides by 255 so that RGB values are in the range 0-1 and the则此函数将其转换为double并除以255，以使RGB值在0-1范围内，
% transformed image can have the proper negative values for a and b.  并且转换后的图像可以具有a和b的正确负值。
%
% See also: LAB2RGB, RGB2NRGB, RGB2CMYK

% Peter Kovesi
% Centre for Exploration Targeting
% The University of Western Australia
% peter.kovesi at uwa edu au

% PK May 2009

function Lab = rgb2lab(im, wp)

    if ~exist('wp', 'var'), wp = 'D65'; end
    
    cform = makecform('srgb2lab',...
                      'adaptedwhitepoint', whitepoint(wp));    
    if strcmp(class(im),'uint8')
        im = double(im)/255;
    end
    Lab = applycform(im, cform);