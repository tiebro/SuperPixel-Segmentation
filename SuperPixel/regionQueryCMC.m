% Find indices of all superpixels adjacent to superpixel n with mean colour
% difference less than Ec.  Use CMC colour difference measure
%   找到与超像素n相邻的所有超像素的索引，其平均色差小于Ec。 使用CMC色差测量
% Arguments:
%             Sp - The struct array of superpixel attributes 超像素属性的struct数组
%             An - Adjacency matrix 邻接矩阵
%              n - Index of superpixel being considered 正在考虑超像素的索引
%             Ec - Colour distance threshold 色彩距离阈值

function neighbours = regionQueryCMC(Sp, Am, n, Ec)
    
    lw = 1;
    neighbours = [];
    
    % Get indices of all superpixels connected to superpixel n  获取连接到超像素n的所有超像素的索引
    ind = find(Am(n,:));
    
    for i = ind

        dE = cmcdifference([Sp(i).L; Sp(i).a; Sp(i).b],...
                           [Sp(n).L; Sp(n).a; Sp(n).b], lw);
        
        if dE < Ec
            neighbours = [neighbours i];     
        end
    end