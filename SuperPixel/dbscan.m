% DBSCAN DBSCAN clustering algorithm
%
% Usage:  [C, ptsC, centres] = dbscan(P, E, minPts)
%
% Arguments:
%         P - dim x Npts array of points.维数x Npts点数组。
%         E - Distance threshold.距离阈值。
%    minPts - Minimum number of points required to form a cluster.形成群集所需的最小点数。
%
% Returns:
%         C - Cell array of length Nc listing indices of points associated with
%             each cluster.    长度为Nc的单元阵列，列出与每个簇相关联的点的索引。
%      ptsC - Array of length Npts listing the cluster number associated with列出与每个点关联的簇编号的长度
%             each point.  If a point is denoted as noise (not enough nearby为Npts的数组。 如果一个点表示
%             elements to form a cluster) its cluster number is 0.为噪声（附近的元素不足以形成一个簇），则其簇号为0。
%   centres - dim x Nc array of the average centre of each cluster.   dim x每个簇的平均中心的Nc数组。

% Reference:
% Martin Ester, Hans-Peter Kriegel, J鰎g Sander, Xiaowei Xu (1996). "A
% density-based algorithm for discovering clusters in large spatial databases
% with noise".  Proceedings of the Second International Conference on Knowledge
% Discovery and Data Mining (KDD-96). AAAI Press. pp. 226-231.  
% Also see: http://en.wikipedia.org/wiki/DBSCAN

% Copyright (c) 2013 Peter Kovesi
% Centre for Exploration Targeting
% The University of Western Australia
% peter.kovesi at uwa edu au
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in 
% all copies or substantial portions of the Software.
%
% The Software is provided "as is", without warranty of any kind.

% PK January 2013

function [C, ptsC, centres] = dbscan(P, E, minPts)
    
    [dim, Npts] = size(P);
    
    ptsC  = zeros(Npts,1);
    C     = {};
    Nc    = 0;               % Cluster counter.
    Pvisit = zeros(Npts,1);  % Array to keep track of points that have been visited.用于跟踪已访问点的数组。
    
    for n = 1:Npts
       if ~Pvisit(n)                            % If this point not visited yet
           Pvisit(n) = 1;                       % mark as visited标记为已访问
           neighbourPts = regionQuery(P, n, E); % and find its neighbours

           if length(neighbourPts) < minPts-1  % Not enough points to form a cluster没有足够的点来形成集群
               ptsC(n) = 0;                    % Mark point n as noise.将点n标记为噪声
           
           else                % Form a cluster...
               Nc = Nc + 1;    % Increment number of clusters and process
                               % neighbourhood.增加集群数和进程邻域数
           
               C{Nc} = [n];    % Initialise cluster Nc with point n使用点n初始化集群Nc
               ptsC(n) = Nc;   % and mark point n as being a member of cluster Nc.并将点n标记为集群Nc的成员。
               
               ind = 1;        % Initialise index into neighbourPts array.将索引初始化为邻居Pts数组。
               
               % For each point P' in neighbourPts ...
               while ind <= length(neighbourPts)
                   
                   nb = neighbourPts(ind);
                   
                   if ~Pvisit(nb)        % If this neighbour has not been visited
                       Pvisit(nb) = 1;   % mark it as visited.
                       
                       % Find the neighbours of this neighbour and if it has找到该邻居的邻居，如果有足够的
                       % enough neighbours add them to the neighbourPts list邻居将其添加到neighbourPts列表中
                       neighbourPtsP = regionQuery(P, nb, E);
                       if length(neighbourPtsP) >= minPts
                           neighbourPts = [neighbourPts  neighbourPtsP];
                       end
                   end            
                   
                   % If this neighbour nb not yet a member of any cluster add it
                   % to this cluster.   如果此邻居nb尚未成为任何群集的成员，则将其添加到此群集
                   if ~ptsC(nb)  
                       C{Nc} = [C{Nc} nb];
                       ptsC(nb) = Nc;
                   end
                   
                   ind = ind + 1;  % Increment neighbour point index and process
                                   % next neighbour 增加邻居点索引并处理下一个邻居
               end
           end
       end
    end
    
    % Find centres of each cluster
    centres = zeros(dim,length(C));
    for n = 1:length(C)
        for k = 1:length(C{n})
            centres(:,n) = centres(:,n) + P(:,C{n}(k));
        end
        centres(:,n) = centres(:,n)/length(C{n});
    end

end % of dbscan    
    
%------------------------------------------------------------------------
% Find indices of all points within distance E of point with index n     查找具有索引n的点的距离E内所有点的索引
% This function could make use of a precomputed distance table to avoid  此功能可以利用预先计算的距离表来避免
% repeated distance calculations, however this would require N^2 storage.  重复距离计算，但是这将需要N ^ 2的存储。
% Not a big problem either way if the number of points being clustered is  如果聚集的点数是
% small.   For large datasets this function will need to be optimised.    对于大型数据集，此功能将需要优化。

% Arguments:
%              P - the dim x Npts array of data points 数据点的dim x Npts数组
%              n - Index of point of interest  兴趣点索引
%              E - Distance threshold  距离阈值

function neighbours = regionQuery(P, n, E)
    
    E2 = E^2;   
    [dim, Npts] = size(P);
    neighbours = [];
    
    for i = 1:Npts
        if i ~= n
            % Test if distance^2 < E^2 
            v = P(:,i)-P(:,n);
            dist2 = v'*v;
            if dist2 < E2 
               neighbours = [neighbours i];     
            end
        end
    end
    
end % of regionQuery

