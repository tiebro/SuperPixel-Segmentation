%--------------------------------------------------------------------------
% DIGIPTS - digitise points in an image���ֻ�ͼ���еĵ�
%
% Function to digitise points in an image.  Points are digitised by clicking�������ֻ�ͼ���еĵ�Ĺ��ܡ�  
% with the left mouse button.  Clicking any other button terminates theͨ�������������������ֻ���
% function.  Each location digitised is marked with a red '+'.�����κ�������ť����ֹ�ù��ܡ� ���ֻ���ÿ��λ�ö����к�ɫ��+����
%
% Usage:  [u,v] = digipts
%
% where u and v are  nx1 arrays of x and y coordinate values digitised in
% the image.  ����u��v��ͼ�������ֻ���x��y����ֵ��nx1�����顣
%
% This function uses the cross-hair cursor provided by GINPUT.  This is
% much more useable than IMPIXEL �˺���ʹ��GINPUT�ṩ��ʮ�ֹ�ꡣ ���IMPIXEL������

function [u,v] = digipts
    
    hold on
    u = []; v = [];
    but = 1;
    while but == 1
	[x y but] = ginput(1);
	if but == 1
	    u = [u;x];
	    v = [v;y];
	    
	    plot(u,v,'r+');
	end
    end
    
    hold off

