%------------------------------------------------------------------------
% Find indices of all superpixels adjacent to superpixel n with mean 
% Euclidean colour difference less than Ec.找到与超像素n相邻的所有超像素的指数，其中平均欧几里德色差小于Ec。
%
% Arguments:
%             Sp - The struct array of superpixel attributes超像素属性的struct数组
%             An - Adjacency matrix邻接矩阵
%              n - Index of point of interest正在考虑超像素的索引
%             Ec - Colour distance threshold色彩距离阈值

function neighbours = regionQueryE(Sp, Am, n, Ec)
    
    E2 = Ec^2;   
    neighbours = [];
    
    % Get indices of all superpixels connected to superpixel n  获取连接到超像素n的所有超像素的索引
    ind = find(Am(n,:));
    
    for i = ind
        % Test if distance^2 < E^2 
        v = Sp(i).value - Sp(n).value;

        if v'*v < E2 
            neighbours = [neighbours i];     
        end
    end
    
    

