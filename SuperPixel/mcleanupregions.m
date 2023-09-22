% MCLEANUPREGIONS  Morphological clean up of small segments in an image of segmented regions
%                  在分割区域的图像中形成小段的形态清理
% Usage: [seg, Am] = mcleanupregions(seg, seRadius)
%
% Arguments: seg - A region segmented image, such as might be produced by a  
%                  graph cut algorithm.  All pixels in each region are labeled
%                  by an integer.  区域分割图像，例如可能由图切割算法产生。每个区域中的所有像素都用整数标记。
%       seRadius - Structuring element radius.  This can be set to 0 in which 结构元素半径。这可以设置为0，
%                  case  the function will simply ensure all labeled regions  在这种情况下，函数将简单地确保
%                  are distinct and relabel them if necessary. 所有标记区域是不同的，并在必要时重新标记它们。
%
% Returns:   seg - The updated segment image.更新的分割图像
%             Am - Adjacency matrix of segments.  Am(i, j) indicates whether 分割的邻接矩阵。 Am（i，j）
%                  segments labeled i and j are connected/adjacent     表示标记为i和j的段是否连接/相邻
%
% Typical application:典型应用
% If a graph cut or superpixel algorithm fails to converge stray segments如果图形切割或超像素算法无法收敛，则可以
% can be left in the result.  This function tries to clean things up by:在结果中保留杂散段。此函数尝试通过以下方式清理：
% 1) Checking there is only one region for each segment label. If there is检查每个段标签只有一个区域。
%    more than one region they are given unique labels. 如果有多个区域，则会为它们提供唯一标签。
% 2) Eliminating regions below the structuring element size 消除结构元素大小以下的区域
%
% Note that regions labeled 0 are treated as a 'privileged' background region请注意，标记为0
% and is not processed/affected by the function.的区域被视为“特权”背景区域，并且不受该函数的处理/影响。
%
% See also: REGIONADJACENCY, RENUMBERREGIONS, CLEANUPREGIONS, MAKEREGIONSDISTINCT

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
%
% March   2013 
% June    2013  Improved morphological cleanup process using distance map

function [seg, Am, mask] = mcleanupregions(seg, seRadius)
option = 2;
    % 1) Ensure every segment is distinct 
    [seg, maxlabel] = makeregionsdistinct(seg);%区分不同区域后的，maxlabel最大标签号
    
    % 2) Perform a morphological opening on each segment, subtract the opening在每个线段
    % from the orignal segment to obtain regions to be reassigned to 上执行形态开口，
    % neighbouring segments.  从原始线段中减去开口，以获得要分配给相邻线段的区域。
    if seRadius
        se = circularstruct(seRadius);   % Accurate and not noticeably slower
                                         % if radius is small 如果半径很小，则准确且不会明显变慢
                                         %得到一个3*3的逻辑矩阵
%       se = strel('disk', seRadius, 4);  % Use approximated disk for speed使用近似磁盘来提高速度
        mask = zeros(size(seg));

        if option == 1        
            for l = 1:maxlabel
                b = seg == l;
                mask = mask | (b - imopen(b,se));%开操作一般使对象的轮廓变得光滑，断开狭窄的间断和消除细的突出物 。
            end
            
        else   % Rather than perform a morphological opening on every 不是按顺序在每个单独的区域上执行
               % individual region in sequence the following finds separate 形态开口，而是在下面找到单独
               % lists of unconnected regions and performs openings on these. 的未连接区域列表并在这些区域
               % Typically an image can be covered with only 5 or 6 lists of 上执行开口。通常，图像可以仅用
               % unconnected regions.  Seems to be about 2X speed of option 5或6个未连接区域列表覆盖。 
               % 1. (I was hoping for more...) 似乎是选项1的2倍速度。（我希望更多......）
            list = finddisconnected(seg);
            
            for n = 1:length(list)
                b = zeros(size(seg));
                for m = 1:length(list{n})
                    b = b | seg == list{n}(m);
                end

                mask = mask | (b - imopen(b,se));
            end
        end
        
        % Compute distance map on inverse of mask计算掩模反转的距离图
        [~, idx] = bwdist(~mask);%bwdist函数用于计算元素之间的距离。
                                %idx表示在该元素所靠近的最近的非零元的位置（即索引值）
        
        % Assign a label to every pixel in the masked area using the label of使用bwdist计算的不在掩码中的最
        % the closest pixel not in the mask as computed by bwdist近像素的标签，为掩模区域中的每个像素分配标签
        seg(mask) = seg(idx(mask));
    end
    
    % 3) As some regions will have been relabled, possibly broken into several由于一些区域将被重新分解，可能被分成
    % parts, or absorbed into others and no longer exist we ensure all regions几个部分，或被吸收到其他区域并且不再存在，
    % are distinct again, and renumber the regions so that they sequentially我们确保所有区域再次不同，并重新编号区域，
    % increase from 1.  We also need to reconstruct the adjacency matrix to使它们从1开始依次增加。我们还需要重建
    % reflect the changed number of regions and their relabeling.邻接矩阵 反映地区数量的变化及其重新标记。

    seg = makeregionsdistinct(seg);
    [seg, minLabel, maxLabel] = renumberregions(seg);
    Am = regionadjacency(seg);    
    
