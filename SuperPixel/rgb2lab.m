% RGB2LAB - RGB to L*a*b* colour space
%
% Usage: Lab = rgb2lab(im, wp)
%
% Arguments:  im - RGB image or Nx3 colourmap for conversion����ת����RGBͼ���Nx3ɫ��ͼ
%             wp - Optional string specifying the adapted white point.��ѡ�ַ�����ָ����Ӧ�İ׵㡣
%                  This defaults to 'D65'.Ĭ��Ϊ��D65��
%
% Returns:   Lab - The converted image or colourmap.ת�����ͼ�����ɫͼ
%
% This function wraps up calls to MAKECFORM and APPLYCFORM in a convenient�˺����Է������ʽ������MAKECFORM��APPLYCFORM�ĵ��á� 
% form.  Note that if the image is of type uint8 this function casts it to��ע�⣬���ͼ�������Ϊuint8��
% double and divides by 255 so that RGB values are in the range 0-1 and the��˺�������ת��Ϊdouble������255����ʹRGBֵ��0-1��Χ�ڣ�
% transformed image can have the proper negative values for a and b.  ����ת�����ͼ����Ծ���a��b����ȷ��ֵ��
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