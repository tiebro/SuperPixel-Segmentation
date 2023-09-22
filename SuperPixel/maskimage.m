% MASKIMAGE Apply mask to image  ��ģͼ��
%
% Usage: maskedim = maskimage(im, mask, col)
%
% Arguments:    im  - Image to be masked Ҫ�ڸǵ�ͼ��
%             mask  - Binary masking image�������ڱ�ͼ��
%              col  - Value/colour to be applied to regions where mask == 1ҪӦ����mask== 1�������ֵ/
%                     If im is a colour image col can be a 3-vector��ɫ���im�ǲ�ɫͼ��
%                     specifying the colour values to be applied.��col������ָ��ҪӦ�õ���ɫֵ��3������
%
% Returns: maskedim - The masked image
%
% See also; DRAWREGIONBOUNDARIES

% Peter Kovesi
% www.peterkovesi.com/matlabfns/
%
% Feb 2013

function maskedim = maskimage(im, mask, col)
    
    [rows,cols, chan] = size(im);
    
    % Set default colour to 0 (black)��Ĭ����ɫ����Ϊ0����ɫ��
    if ~exist('col', 'var'), col = 0; end
    
    % Ensure col has same length as image depth.ȷ��col�ĳ�����ͼ�������ͬ��
    if length(col) == 1
        col = repmat(col, [chan 1]);
    else
        assert(length(col) == chan);
    end
    
    % Perform maskingִ����ģ
    maskedim = im;
    for n = 1:chan
        tmp = maskedim(:,:,n);
        tmp(mask) = col(n);
        maskedim(:,:,n) = tmp;
    end
    