% MASKIMAGE Apply mask to image  掩模图像
%
% Usage: maskedim = maskimage(im, mask, col)
%
% Arguments:    im  - Image to be masked 要掩盖的图像
%             mask  - Binary masking image二进制掩蔽图像
%              col  - Value/colour to be applied to regions where mask == 1要应用于mask== 1的区域的值/
%                     If im is a colour image col can be a 3-vector颜色如果im是彩色图像，
%                     specifying the colour values to be applied.则col可以是指定要应用的颜色值的3向量。
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
    
    % Set default colour to 0 (black)将默认颜色设置为0（黑色）
    if ~exist('col', 'var'), col = 0; end
    
    % Ensure col has same length as image depth.确保col的长度与图像深度相同。
    if length(col) == 1
        col = repmat(col, [chan 1]);
    else
        assert(length(col) == chan);
    end
    
    % Perform masking执行掩模
    maskedim = im;
    for n = 1:chan
        tmp = maskedim(:,:,n);
        tmp(mask) = col(n);
        maskedim(:,:,n) = tmp;
    end
    