% TESTDBSCAN   Program to test/demonstrate the DBSCAN clustering algorithm���ڲ���/��ʾDBSCAN�����㷨�ĳ���
%
% Simple usage:             testdbscan;
%
% Full usage:   [C, ptsC] = testdbscan(E, minPts)
%       
% 
% Arguments:    
%         E - Distance threshold for clustering. Defaults to 0.3    E-Ⱥ���ľ�����ֵ��Ĭ��Ϊ0.3
%    minPts - Minimum number of points required to form a cluster. 
%             Defaults to 3    �γ�Ⱥ���������С������ Ĭ��Ϊ3
%
% Returns:
%         C - Cell array of length Nc listing indices of points associated with
%             each cluster.      ����ΪNc�ĵ�Ԫ���У��г���ÿ����������ĵ��������
%      ptsC - Array of length Npts listing the cluster number associated with�г���ÿ��������Ĵر�ŵ�
%             each point.  If a point is denoted as noise (not enough nearby����ΪNpts�����顣 ���һ�����ʾΪ
%             elements to form a cluster) its cluster number is 0.  ������������Ԫ�ز������γ�һ���أ�������غ�Ϊ0��
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
    
    % Perform clusteringִ��Ⱥ��
    P = [x'; y'];
    [C, ptsC, centres] = dbscan(P, E, minPts);
    
    
    for n = 1:length(x)
        text(x(n),y(n)+.04, sprintf('%d',ptsC(n)), 'color', [0 0 1]);
    end
    title('Points annotated by cluster number')
    hold off
    
% %--------------------------------------------------------------------------
% % DIGIPTS - digitise points in an image���ֻ�ͼ���еĵ�
% %
% % Function to digitise points in an image.  Points are digitised by clicking�������ֻ�ͼ���еĵ�Ĺ��ܡ�  
% % with the left mouse button.  Clicking any other button terminates theͨ�������������������ֻ���
% % function.  Each location digitised is marked with a red '+'.�����κ�������ť����ֹ�ù��ܡ� ���ֻ���ÿ��λ�ö����к�ɫ��+����
% %
% % Usage:  [u,v] = digipts
% %
% % where u and v are  nx1 arrays of x and y coordinate values digitised in
% % the image.  ����u��v��ͼ�������ֻ���x��y����ֵ��nx1�����顣
% %
% % This function uses the cross-hair cursor provided by GINPUT.  This is
% % much more useable than IMPIXEL �˺���ʹ��GINPUT�ṩ��ʮ�ֹ�ꡣ ���IMPIXEL������
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
