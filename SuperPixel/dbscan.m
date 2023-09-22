% DBSCAN DBSCAN clustering algorithm
%
% Usage:  [C, ptsC, centres] = dbscan(P, E, minPts)
%
% Arguments:
%         P - dim x Npts array of points.ά��x Npts�����顣
%         E - Distance threshold.������ֵ��
%    minPts - Minimum number of points required to form a cluster.�γ�Ⱥ���������С������
%
% Returns:
%         C - Cell array of length Nc listing indices of points associated with
%             each cluster.    ����ΪNc�ĵ�Ԫ���У��г���ÿ����������ĵ��������
%      ptsC - Array of length Npts listing the cluster number associated with�г���ÿ��������Ĵر�ŵĳ���
%             each point.  If a point is denoted as noise (not enough nearbyΪNpts�����顣 ���һ�����ʾ
%             elements to form a cluster) its cluster number is 0.Ϊ������������Ԫ�ز������γ�һ���أ�������غ�Ϊ0��
%   centres - dim x Nc array of the average centre of each cluster.   dim xÿ���ص�ƽ�����ĵ�Nc���顣

% Reference:
% Martin Ester, Hans-Peter Kriegel, J�rg Sander, Xiaowei Xu (1996). "A
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
    Pvisit = zeros(Npts,1);  % Array to keep track of points that have been visited.���ڸ����ѷ��ʵ�����顣
    
    for n = 1:Npts
       if ~Pvisit(n)                            % If this point not visited yet
           Pvisit(n) = 1;                       % mark as visited���Ϊ�ѷ���
           neighbourPts = regionQuery(P, n, E); % and find its neighbours

           if length(neighbourPts) < minPts-1  % Not enough points to form a clusterû���㹻�ĵ����γɼ�Ⱥ
               ptsC(n) = 0;                    % Mark point n as noise.����n���Ϊ����
           
           else                % Form a cluster...
               Nc = Nc + 1;    % Increment number of clusters and process
                               % neighbourhood.���Ӽ�Ⱥ���ͽ���������
           
               C{Nc} = [n];    % Initialise cluster Nc with point nʹ�õ�n��ʼ����ȺNc
               ptsC(n) = Nc;   % and mark point n as being a member of cluster Nc.������n���Ϊ��ȺNc�ĳ�Ա��
               
               ind = 1;        % Initialise index into neighbourPts array.��������ʼ��Ϊ�ھ�Pts���顣
               
               % For each point P' in neighbourPts ...
               while ind <= length(neighbourPts)
                   
                   nb = neighbourPts(ind);
                   
                   if ~Pvisit(nb)        % If this neighbour has not been visited
                       Pvisit(nb) = 1;   % mark it as visited.
                       
                       % Find the neighbours of this neighbour and if it has�ҵ����ھӵ��ھӣ�������㹻��
                       % enough neighbours add them to the neighbourPts list�ھӽ�����ӵ�neighbourPts�б���
                       neighbourPtsP = regionQuery(P, nb, E);
                       if length(neighbourPtsP) >= minPts
                           neighbourPts = [neighbourPts  neighbourPtsP];
                       end
                   end            
                   
                   % If this neighbour nb not yet a member of any cluster add it
                   % to this cluster.   ������ھ�nb��δ��Ϊ�κ�Ⱥ���ĳ�Ա��������ӵ���Ⱥ��
                   if ~ptsC(nb)  
                       C{Nc} = [C{Nc} nb];
                       ptsC(nb) = Nc;
                   end
                   
                   ind = ind + 1;  % Increment neighbour point index and process
                                   % next neighbour �����ھӵ�������������һ���ھ�
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
% Find indices of all points within distance E of point with index n     ���Ҿ�������n�ĵ�ľ���E�����е������
% This function could make use of a precomputed distance table to avoid  �˹��ܿ�������Ԥ�ȼ���ľ����������
% repeated distance calculations, however this would require N^2 storage.  �ظ�������㣬�����⽫��ҪN ^ 2�Ĵ洢��
% Not a big problem either way if the number of points being clustered is  ����ۼ��ĵ�����
% small.   For large datasets this function will need to be optimised.    ���ڴ������ݼ����˹��ܽ���Ҫ�Ż���

% Arguments:
%              P - the dim x Npts array of data points ���ݵ��dim x Npts����
%              n - Index of point of interest  ��Ȥ������
%              E - Distance threshold  ������ֵ

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

