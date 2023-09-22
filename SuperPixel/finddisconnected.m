% FINDDISCONNECTED find groupings of disconnected labeled regions
%                   找到断开连接的标记区域的分组
% Usage: list = finddisconnected(l)
%
% Argument:   l - A labeled image segmenting an image into regions, such as  标记图像将图像分割成区域，例如
%                 might be produced by a graph cut or superpixel algorithm. 可以通过图形切割或超像素算法产生。
%                 All pixels in each region are labeled by an integer.      每个区域中的所有像素都用整数标记。
%
% Returns: list - A cell array of lists of regions that are not未连接的区域列表的单元阵列。 
%                 connected. Typically there are 5 to 6 lists.通常有5到6个列表。
%
% Used by MCLEANUPREGIONS to reduce the number of morphological closing
% operations  由MCLEANUPREGIONS用于减少形态闭合操作的数量
%
% See also: MCLEANUPREGIONS, REGIONADJACENCY

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

% PK July 2013


function list = finddisconnected(l)
 
    debug = 0;
    [Am, Al] = regionadjacency(l);%计算标记的分割区域的图像的邻接矩阵
    
    N = max(l(:));  % number of labels
    
    % Array for keeping track of visited labels用于跟踪访问过的标签的数组
    visited = zeros(N,1);

    list = {};
    listNo = 0;
    for n = 1:N

        if ~visited(n)
            listNo = listNo + 1;
            list{listNo} = n;
            visited(n) = 1;
            
            % Find all regions not directly connected to n and not visited查找未直接连接到n且未访问的所有区域
            notConnected = setdiff(find(~Am(n,:)), find(visited));
            %c = setdiff(A, B)   返回在A中有，而B中没有的值，结果向量将以升序排序返回。
            
            % For each unconnected region check that it is not already对于每个未连接的区域，检查它是否尚未
            % connected to a region in the list. If not, add to list连接到列表中的区域。 如果没有，请添加到列表中
            for m = notConnected
                if isempty(intersect(Al{m}, list{listNo}))
                    list{listNo} = [list{listNo} m];
                    visited(m) = 1;
                end
            end
         end % if not visited(n)
        
    end
    
    % Display each list of unconncted regions as an image 将每个未连接区域列表显示为图像
    if debug   %调试模式
        for n = 1:length(list)
            
            mask = zeros(size(l));
            for m = 1:length(list{n})
                mask = mask | l == list{n}(m);
            end
            
            fprintf('list %d of %d length %d \n', n, length(list), length(list{n}))
            show(mask);
            keypause
        end
    end
    