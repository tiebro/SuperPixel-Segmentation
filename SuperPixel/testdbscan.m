% TESTDBSCAN   Program to test/demonstrate the DBSCAN clustering algorithm用于测试/演示DBSCAN聚类算法的程序
%
% Simple usage:             testdbscan;
%
% Full usage:   [C, ptsC] = testdbscan(E, minPts)
%       
% 
% Arguments:    
%         E - Distance threshold for clustering. Defaults to 0.3    E-群集的距离阈值。默认为0.3
%    minPts - Minimum number of points required to form a cluster. 
%             Defaults to 3    形成群集所需的最小点数。 默认为3
%
% Returns:
%         C - Cell array of length Nc listing indices of points associated with
%             each cluster.      长度为Nc的单元阵列，列出与每个簇相关联的点的索引。
%      ptsC - Array of length Npts listing the cluster number associated with列出与每个点关联的簇编号的
%             each point.  If a point is denoted as noise (not enough nearby长度为Npts的数组。 如果一个点表示为
%             elements to form a cluster) its cluster number is 0.  噪声（附近的元素不足以形成一个簇），则其簇号为0。
%
% See also: DBSCAN

% Copyright (c) 2013 Peter Kovesi
% Centre for Exploration Targeting
% The University of Western Australia
% peter.kovesi at uwa edu au
%
% Jan 2013

function [C, ptsC, centres] = testdbscan(E, minPts)
    
    if ~exist('E', 'var'), E = 0.3; end;
    if ~exist('minPts', 'var'), minPts = 3; end;
    
    figure(1), clf, axis([-1 1 -1 1]);
    
    fprintf('Digitise a series of points that form some clusters.  Right-click to finish\n');
    [x,y] = digipts;
    hold on
    
    % Perform clustering执行群集
    P = [x'; y'];
    [C, ptsC, centres] = dbscan(P, E, minPts);
    
    
    for n = 1:length(x)
        text(x(n),y(n)+.04, sprintf('%d',ptsC(n)), 'color', [0 0 1]);
    end
    title('Points annotated by cluster number')
    hold off
    
% %--------------------------------------------------------------------------
% % DIGIPTS - digitise points in an image数字化图像中的点
% %
% % Function to digitise points in an image.  Points are digitised by clicking用于数字化图像中的点的功能。  
% % with the left mouse button.  Clicking any other button terminates the通过单击鼠标左键将点数字化。
% % function.  Each location digitised is marked with a red '+'.单击任何其他按钮将终止该功能。 数字化的每个位置都标有红色“+”。
% %
% % Usage:  [u,v] = digipts
% %
% % where u and v are  nx1 arrays of x and y coordinate values digitised in
% % the image.  其中u和v是图像中数字化的x和y坐标值的nx1个数组。
% %
% % This function uses the cross-hair cursor provided by GINPUT.  This is
% % much more useable than IMPIXEL 此函数使用GINPUT提供的十字光标。 这比IMPIXEL更有用
% 
% function [u,v] = digipts
%     
%     hold on
%     u = []; v = [];
%     but = 1;
%     while but == 1
% 	[x y but] = ginput(1);
% 	if but == 1
% 	    u = [u;x];
% 	    v = [v;y];
% 	    
% 	    plot(u,v,'r+');
% 	end
%     end
%     
%     hold off
