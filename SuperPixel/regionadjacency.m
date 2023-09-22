% REGIONADJACENCY Computes adjacency matrix for image of labeled segmented regions
%                 计算标记的分割区域的图像的邻接矩阵
% Usage:  [Am, Al] = regionadjacency(L, connectivity)
%
% Arguments:  L - A region segmented image, such as might be produced by a  区域分割图像，例如
%                 graph cut or superpixel algorithm.  All pixels in each    可以通过图形切割或超像素
%                 region are labeled by an integer.        算法产生。 每个区域中的所有像素都用整数标记。
%  connectivity - 8 or 4.  If not specified connectivity defaults to 8. 如果未指定，则连接默认为8。
%
% Returns:   Am - An adjacency matrix indicating which labeled regions are邻接矩阵，指示哪些
%                 adjacent to each other, that is, they share boundaries. Am标记区域彼此相邻，
%                 is sparse to save memory.即，它们共享边界。 为了节省内存，Am很稀疏。
%            Al - A cell array representing the adjacency list corresponding
%                 to Am.  Al{n} is an array of the region indices adjacent to
%                 region n.    表示与Am对应的邻接列表的单元阵列。 Al{n}是与区域n相邻的区域索引的数组。
%
% Regions with a label of 0 are not processed. They are considered to be
% 'background regions' that are not to be considered.  If you want to include
% these regions you should assign a new positive label to these areas using, say
% >> L(L==0) = max(L(:)) + 1;            标签为0的区域不会被处理。 它们被认为是不被考虑的“背景区域”。
% 如果你想要包含这些区域，你应该使用，例如>> L（L == 0）= max（L（:)）+ 1为这些区域分配一个新的正面标签;
% See also: CLEANUPREGIONS, RENUMBERREGIONS, SLIC

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
% February 2013  Original version
% July     2013  Speed improvement in sparse matrix formation (4x)

function  [Am, varargout] = regionadjacency(L, connectivity)

    if ~exist('connectivity', 'var'), connectivity = 4; end
    [rows,cols] = size(L);
    
    % Identify the unique labels in the image, excluding 0 as a label.识别图像中的唯一标签，不包括0作为标签。
    labels = setdiff(unique(L(:))',0);%L(:)中不为0的数按大小排列

    if isempty(labels) %判断是否为空
        warning('There are no objects in the image')
        Am = [];
        Al = {};
        return
    end

    N = max(labels);    % Required size of adjacency matrix邻接矩阵的大小
    
    % Strategy:  Step through the labeled image.  For 8-connectedness inspect 策略：逐步完成标记图像。 用于8连通性检查
    % pixels as follows and set the appropriate entries in the adjacency
    % matrix. 像素如下，并在邻接矩阵中设置适当的条目。
    %      x - o
    %    / | \
    %  o   o   o
    %
    % For 4-connectedness we only inspect the following pixels对于4连通性，我们只检查以下像素
    %      x - o
    %      | 
    %      o  
    %
    % Becuase the adjacency search looks 'forwards' a final OR operation is 因为邻接搜索
    % performed on the adjacency matrix and its transpose to ensure  看起来是“向前”，
    % connectivity both ways.  所以对邻接矩阵及其转置执行最终OR运算以确保两种方式的连接。

    % Allocate vectors for forming row, col, value triplets used to construct分配用于形成用于构造
    % sparse matrix.  Forming these vectors first is faster than filling稀疏矩阵的行，列，值三元组
    % entries directly into the sparse matrix的向量。首先形成这些向量比将条目直接填充到稀疏矩阵中更快
    i = zeros(rows*cols,1);  % row value
    j = zeros(rows*cols,1);  % col value
    s = zeros(rows*cols,1);  % value
    
    if connectivity == 8
        n = 1;
        for r = 1:rows-1

            % Handle pixels in 1st column处理第1列中的像素
            i(n) = L(r,1); j(n) = L(r  ,2); s(n) = 1; n=n+1;
            i(n) = L(r,1); j(n) = L(r+1,1); s(n) = 1; n=n+1;
            i(n) = L(r,1); j(n) = L(r+1,2); s(n) = 1; n=n+1;
            
            % ... now the rest of the column
            for c = 2:cols-1
               i(n) = L(r,c); j(n) = L(r  ,c+1); s(n) = 1; n=n+1;
               i(n) = L(r,c); j(n) = L(r+1,c-1); s(n) = 1; n=n+1;
               i(n) = L(r,c); j(n) = L(r+1,c  ); s(n) = 1; n=n+1;
               i(n) = L(r,c); j(n) = L(r+1,c+1); s(n) = 1; n=n+1;
            end
        end
        
    elseif connectivity == 4
        n = 1;
        for r = 1:rows-1
            for c = 1:cols-1
                i(n) = L(r,c); j(n) = L(r  ,c+1); s(n) = 1; n=n+1;
                i(n) = L(r,c); j(n) = L(r+1,c  ); s(n) = 1; n=n+1;
            end
        end
    
    else
        error('Connectivity must be 4 or 8');
    end
    
    % Form the logical sparse adjacency matrix形成逻辑稀疏邻接矩阵
    Am = logical(sparse(i, j, s, N, N)); 
    
    % Zero out the diagonal 将对角线清零
    for r = 1:N
        Am(r,r) = 0;
    end
    
    % Ensure connectivity both ways for all regions.确保所有地区的双向连接。
    Am = Am | Am';
    
    % If an adjacency list is requested...如果要求邻接列表......
    if nargout == 2
        Al = cell(N,1);
        for r = 1:N
            Al{r} = find(Am(r,:));
        end
        varargout{1} = Al;
    end
    
