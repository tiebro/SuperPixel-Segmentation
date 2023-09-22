%--------------------------------------------------------------------------
% DIGIPTS - digitise points in an image数字化图像中的点
%
% Function to digitise points in an image.  Points are digitised by clicking用于数字化图像中的点的功能。  
% with the left mouse button.  Clicking any other button terminates the通过单击鼠标左键将点数字化。
% function.  Each location digitised is marked with a red '+'.单击任何其他按钮将终止该功能。 数字化的每个位置都标有红色“+”。
%
% Usage:  [u,v] = digipts
%
% where u and v are  nx1 arrays of x and y coordinate values digitised in
% the image.  其中u和v是图像中数字化的x和y坐标值的nx1个数组。
%
% This function uses the cross-hair cursor provided by GINPUT.  This is
% much more useable than IMPIXEL 此函数使用GINPUT提供的十字光标。 这比IMPIXEL更有用

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

