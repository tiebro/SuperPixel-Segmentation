% MAKEREGIONSDISTINCT Ensures labeled segments are distinct确保标记的段是不同的
%
% Usage: [seg, maxlabel] = makeregionsdistinct(seg, connectivity)
%
% Arguments: seg - A region segmented image, such as might be produced by a  区域分割图像， 
%                  superpixel or graph cut algorithm.  All pixels in each    例如可能由超像素或图形切割
%                  region are labeled by an integer.           算法产生。每个区域中的所有像素都用整数标记。
%   connectivity - Optional parameter indicating whether 4 or 8 connectedness
%                  should be used.  Defaults to 4. 可选参数，指示是否应使用4或8连通性。 默认为4。
%
% Returns:   seg - A labeled image where all segments are distinct.标记图像，其中所有段都是不同的。
%       maxlabel - Maximum segment label number.最大段标签号。
%
% Typical application: A graphcut or superpixel algorithm may terminate in a few  图形切割或超像素算法可以
% cases with multiple regions that have the same label.  This function  在具有相同标签的多个区域的少数情况下 
% identifies these regions and assigns a unique label to them.  终止。此功能可识别这些区域并为其分配唯一标签。
%
% See also: SLIC, CLEANUPREGIONS, RENUMBERREGIONS

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

% June 2013


function [seg, maxlabel] = makeregionsdistinct(seg, connectivity)
    
    if ~exist('connectivity', 'var'), connectivity = 4; end     %不存在=4
    
    % Ensure every segment is distinct but do not touch segments 
    % with a label of 0 确保每个段都是不同的，但不要触摸标签为0的段
    labels = unique(seg(:))';%unique取集合的不重复元素构成的向量
    maxlabel = max(labels);
    labels = setdiff(labels,0);  % Remove 0 from the label list从标签列表中删除0
    
    for l = labels
        [bl,num] = bwlabel(seg==l, connectivity);  %L返回seg==l中每个连通区域的类别标签，num返回连通区域的个数
        
        if num > 1  % We have more than one region with the same label我们有多个区域具有相同的标签
            for n = 2:num
                maxlabel = maxlabel+1;  % Generate a new label生成新标签
                seg(bl==n) = maxlabel;  % and assign to this segment并分配给此细分
            end
        end
    end

