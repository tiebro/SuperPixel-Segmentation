% DRAWREGIONBOUNDARIES Draw boundaries of labeled regions in an image     DRAWREGIONBOUNDARIES����ͼ���б������ı߽�
%
% Usage: maskim = drawregionboundaries(l, im, col)
%
% Arguments:����
%            l - Labeled image of regions.�����ǵ�ͼ��
%           im - Optional image to overlay the region boundaries on.��ѡͼ���Ը�������߽�
%          col - Optional colour specification. Defaults to black.  Note that��ѡ��ɫ��� Ĭ��Ϊ��ɫ�� ��ע�⣬
%                the colour components are specified as values 0-255.��ɫ����ָ��Ϊֵ0-255��
%                For example red is [255 0 0] and white is [255 255 255].���磬��ɫ��[255 0 0]����ɫ��[255 255 255]��
%
% Returns: 
%       maskim - If no image has been supplied maskim is a binary mask���û���ṩͼ����maskim�Ƕ�������ģͼ��
%                image indicating where the region boundaries are.ָʾ����߽��λ�á�
%                If an image has been supplied maskim is the image with the����ṩ��ͼ����maskim�Ǹ�������߽��ͼ��
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
    
    % Form the mask by applying a sobel edge detector to the labeled image,ͨ����sobel��Ե�����Ӧ���ڱ��ͼ�����γ���ģ��
    % thresholding and then thinning the result.��ֵ����Ȼ��ϡ������
%    h = [1  0 -1
%         2  0 -2
%         1  0 -1];
    h = [-1 1];  % A simple small filter is better in this application.
                 % Small regions 1 pixel wide get missed using a Sobel
                 % operator һ���򵥵�С�˲��������Ӧ���и��á�ʹ��Sobel���ӿ��Զ�ʧ1�����ؿ��С����
% h = fspecial('prewitt');%�˲�����
    gx = filter2(h ,l);
    gy = filter2(h',l);
    maskim = (gx.^2 + gy.^2) > 0;
    maskim = bwmorph(maskim, 'thin', Inf);%�Զ�ֵͼ��maskim����ϸ������
    
    % Zero out any mask values that may have been set around the edge of the
    % image.  �����������ͼ���Ե��Χ���õ��κ���ģֵ��
    maskim(1,:) = 0; maskim(end,:) = 0;
    maskim(:,1) = 0; maskim(:,end) = 0;
    
    % If an image has been supplied apply the mask to the image and return it ������ṩͼ���뽫��ģӦ����ͼ�񲢽��䷵��
    if exist('im', 'var') 
        if ~exist('col', 'var'), col = 0; end
        maskim = maskimage(im, maskim, col);
    end