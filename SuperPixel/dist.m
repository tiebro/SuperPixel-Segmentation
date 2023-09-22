%-- dist -------------------------------------------
%
% Usage:  D = dist(C, im, r1, c1, S, m)
% 
% Arguments:   C - Cluster being considered正在考虑集群
%             im - sub-image surrounding cluster centre围绕集群中心的子图像
%         r1, c1 - row and column of top left corner of sub image within the
%                  overall image.   整个图像中子图像左上角的行和列。
%              S - grid spacing 网格间距
%              m - weighting factor between colour and spatial differences.颜色和空间差异之间的加权因子。
%
% Returns:     D - Distance image giving distance of every pixel in the
%                  subimage from the cluster centre  距离图像给出子图像中每个像素距离聚类中心的距离
%
% Distance = sqrt( dc^2 + (ds/S)^2*m^2 )
% where:
% dc = sqrt(dl^2 + da^2 + db^2)  % Colour distance
% ds = sqrt(dx^2 + dy^2)         % Spatial distance
%
% m is a weighting factor representing the nominal maximum colour distancem是表示预期的标称最大颜色距离的加权
% expected so that one can rank colour similarity relative to distance因子，因此可以相对于距离相似性对颜色相似性
% similarity.  try m in the range [1-40] for L*a*b* space进行排序。 对于L* a* b*空间，尝试在[1-40]范围内的m
%
% ?? Might be worth trying the Geometric Mean instead ?? 可能值得尝试几何平均值
%  Distance = sqrt(dc * ds)
% but having a factor 'm' to play with is probably handy  但有一个因素'm'可能很方便

% This code could be more efficient  此代码可能更有效

function D = dist(C, im, r1, c1, S, m)

    % Squared spatial distance  平方的空间距离
    %    ds is a fixed 'image' we should be able to exploit this  ds是一个固定的'图像'，我们应该能够
    %    and use a fixed meshgrid for much of the time somehow...  利用它并在某种程度上使用固定的网格...
    [rows, cols, chan] = size(im);
    [x,y] = meshgrid(c1:(c1+cols-1), r1:(r1+rows-1));
    x = x-C(4);  % x and y dist from cluster centre   x和y远离集群中心
    y = y-C(5);
    ds2 = x.^2 + y.^2;
    
    % Squared colour difference平方色差
    for n = 1:3
        im(:,:,n) = (im(:,:,n)-C(n)).^2;
    end
    dc2 = sum(im,3);
    
    D = sqrt(dc2 + ds2/S^2*m^2);

