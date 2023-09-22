% RENUMBERREGIONS
%
% Usage: [nL, minLabel, maxLabel] = renumberregions(L)
%
% Argument:   L - A labeled image segmenting an image into regions, such as标记图像将图像分割成区域，
%                 might be produced by a graph cut or superpixel algorithm.例如可以通过图形切割或超像素算法产生。
%                 All pixels in each region are labeled by an integer.每个区域中的所有像素都用整数标记。
%
% Returns:   nL - A relabeled version of L so that label numbers form aL的重新标记版本，使得标签号
%                 sequence 1:maxRegions  or 0:maxRegions-1 depending on形成序列1：maxRegions或0：maxRegions-1，
%                 whether L has a region labeled with 0s or not.这取决于L是否具有标记为0s的区域。
%      minLabel - Minimum label in the renumbered image.  This will be 0 or 1.重新编号图像中的最小标签。
%      maxLabel - Maximum label in the renumbered image.
%
% Application: Segmentation algorithms can produce a labeled image with a non分割算法可以产生具有区域1 4 6等的
% contiguous numbering of regions 1 4 6 etc. This function renumbers them into a非连续编号的标记图像。该函数
% contiguous sequence.  If the input image has a region labeled with 0s this将它们重新编号为连续序列。 如果输入图像
% region is treated as a privileged 'background region' and retains its 0具有标记为0的区域，则该区域被视为特权“背景区域”
% labeling. The resulting image will have labels ranging over 0:maxRegions-1.并保留其0标记。 生成的图像将具有范围
% Otherwise the image will be relabeled over the sequence 1:maxRegions超过0的标签：maxRegions-1。 否则，图像将在序列1：maxRegions上重新标记
%
% See also: CLEANUPREGIONS, REGIONADJACENCY

% Copyright (c) 2010 Peter Kovesi
% www.peterkovesi.com/matlabfns/
%
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in 
% all copies or substantial portions of the Software.
%
% October  2010
% February 2013 Return label numbering range

function [nL, minLabel, maxLabel] = renumberregions(L)

    nL = L;
    labels = unique(L(:))';  % Sorted list of unique labels已排序的唯一标签列表
    N = length(labels);
    
    % If there is a label of 0 we ensure that we do not renumber that region如果标签为0，我们确保不会
    % by removing it from the list of labels to be renumbered.通过从要重新编号的标签列表中删除该区域来重新编号。
    if labels(1) == 0
        labels = labels(2:end);
        minLabel = 0;
        maxLabel = N-1;
    else
        minLabel = 1;
        maxLabel = N;
    end
    
    % Now do the relabelling现在进行重新标记
    count = 1;
    for n = labels
        nL(L==n) = count;
        count = count+1;
    end
    