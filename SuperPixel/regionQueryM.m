%------------------------------------------------------------------------
% Find indices of all superpixels adjacent to superpixel n with mean Lab 
% colour difference less than Ec. 查找与超像素n相邻的所有超像素的索引，其中Lab色差小于Ec。
%
% Arguments:
%             Sp - The struct array of superpixel attributes超像素属性的struct数组
%             An - Adjacency matrix邻接矩阵
%              n - Index of point of interest正在考虑超像素的索引
%             Ec - Colour distance threshold色彩距离阈值

function neighbours = regionQueryM(Sp, Am, n, Ec, gray, E3)
    
    E2 = Ec^2;
    neighbours = [];
    
    % Get indices of all superpixels connected to superpixel n  获取连接到超像素n的所有超像素的索引
    ind = find(Am(n,:));
    
    for i = ind
        % Test if distance^2 < E^2 
%         v = [Sp(i).L; 1000*Sp(i).a; 100000*Sp(i).b] - ...
%             [Sp(n).L; 1000*Sp(n).a; 100000*Sp(n).b];
        v = [Sp(i).L; Sp(i).a; Sp(i).b] - ...
            [Sp(n).L; Sp(n).a; Sp(n).b];

        dist2 = v'*v;
        
        
%         v1 = [Sp(i).stdL; Sp(i).stda; Sp(i).stdb] - ...
%             [Sp(n).stdL; Sp(n).stda; Sp(n).stdb];
%         dist3 = sum(sum(abs(v1)));
%         dist3 = v1'*v1;     
        
        if dist2 < E2
%         if dist2 < E2 && dist3 < E3
            neighbours = [neighbours i];     
        end


%纹理特征        
%         [row,col] = size(gray);
%         r1 = round(Sp(i).r);
%         c1 = round(Sp(i).c);
%         r2 = round(Sp(n).r);
%         c2 = round(Sp(n).c);
%         MIN1 = min([r1,c1,r2,c2]);
%         MAX1 = max([r1,r2]);
%         MAX2 = max([c1,c2]);
%         if MIN1 > 8 && MAX1 < (row-8) && MAX2 < (col-8)
%         gray1 = gray(r1-8:r1+8,c1-8:c1+8);
%         gray2 = gray(r2-8:r2+8,c2-8:c2+8);
%         [s1, s2, s3, s4, s5] = glcm(gray1);
%         a1 = s1;a2 = s2;a3 = s3;a4 = s4;a5 = s5;
%         [s1, s2, s3, s4, s5] = glcm(gray2);
%         b1 = s1;b2 = s2;b3 = s3;b4 = s4;b5 = s5;
%         
%         d1 = 0;d2 = 0;d3 = 0;d4 = 0;d5 = 0;
%         if abs(a1-b1) < 10
%             d1 = 1;
%         end
%         if abs(a2-b2) < 0.05
%             d2 = 1;
%         end
%         if abs(a3-b3)/abs(a3+b3) < 0.5
%             d3 = 1;
%         end
%         if abs(a4-b4)/abs(a4+b4) < 0.5
%             d4 = 1;
%         end
%         if abs(a5-b5)/abs(a5+b5) < 0.05
%             d5 = 1;
%         end
%         if d1 == 1 && d2 == 1 && d3 == 1 && d4 == 1 && d5 == 1
%             neighbours = [neighbours i];
%         end
%         end
    end

