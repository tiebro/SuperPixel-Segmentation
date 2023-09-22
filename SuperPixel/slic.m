% SLIC Simple Linear Iterative Clustering SuperPixels简单线性迭代聚类超像素
%
% Implementation of Achanta, Shaji, Smith, Lucchi, Fua and Susstrunk's
% SLIC Superpixels
%
% Usage:   [l, Am, Sp, d] = slic(im, k, m, seRadius, colopt, mw)
%
% Arguments:  im - Image to be segmented.要分割的图像
%              k - Number of desired superpixels. Note that this is nominal所需超像素的数量。注意，这是标称的，
%                  the actual number of superpixels generated will generally所生成的超像素的实际数量通常会稍大，
%                  be a bit larger, especially if parameter m is small.特别是如果参数m很小。
%              m - Weighting factor between colour and spatial   颜色和空间差异之间的加权因子。   
%                  differences. Values from about 5 to 40 are useful.  Use a  大约5到40的值是有用的。使用较
%                  large value to enforce superpixels with more regular and大的值来强制使用更规则和更平滑的
%                  smoother shapes. Try a value of 10 to start with.形状的超像素。尝试使用值10开始。
%       seRadius - Regions morphologically smaller than this are merged with  形态学上小于此的区域与相邻区域合并。 
%                  adjacent regions. Try a value of 1 or 1.5.  Use 0 to       尝试值1或1.5。 使用0禁用。
%                  disable.
%         colopt - String 'mean' or 'median' indicating how the cluster  字符串'mean'或'median'表示应如何
%                  colour centre should be computed. Defaults to 'mean'  计算簇颜色中心。 默认为'均值'
%             mw - Optional median filtering window size.  Image compression 可选的中值过滤窗口大小。图像压缩    
%                  can result in noticeable artifacts in the a*b* components 可能会导致图像的a* b*组件出现明显的
%                  of the image.  Median filtering can reduce this. mw can be 瑕疵。中值过滤可以减少这种情况。
%                  a single value in which case the same median filtering ismw 可以是单个值，在这种情况下，
%                  applied to each L* a* and b* components.  Alternatively it 对每个L*a*和b*分量应用相同的中值
%                  can be a 2-vector where mw(1) specifies the median  滤波。或者，它可以是2向量，其中mw（1）
%                  filtering window to be applied to L* and mw(2) is the 指定要应用于L*的中值滤波窗口，
%                  median filtering window to be applied to a* and b*. 而mw（2）是要应用于a*和b*的中值滤波窗口。
%
% Returns:     l - Labeled image of superpixels. Labels range from 1 to k.标签图像的超像素。 标签范围从1到k。
%             Am - Adjacency matrix of segments.  Am(i, j) indicates whether  分割的邻接矩阵。 Am（i，j）
%                  segments labeled i and j are connected/adjacent      表示标记为i和j的段是否连接/相邻
%             Sp - Superpixel attribute structure array with fields:  具有字段的超像素属性结构数组：
%                   L  - Mean L* value 平均L*值
%                   a  - Mean a* value 平均a*值
%                   b  - Mean b* value 平均b*值
%                   r  - Mean row value 平均row值
%                   c  - Mean column value 平均column值
%                   stdL  - Standard deviation of L*     L*的标准偏差
%                   stda  - Standard deviation of a*     a*的标准偏差 
%                   stdb  - Standard deviation of b*     b*的标准偏差
%                   N - Number of pixels
%                   edges - List of edge numbers that bound each  绑定每个超像素的边数列表。 
%                           superpixel. This field is allocated, but not set, 该字段由SLIC分配，
%                           by SLIC. Use  for this.  但未设置。使用SPEDGES。
%              d - Distance image giving the distance each pixel is from its
%                  associated superpixel centre.距离图像给出每个像素距其相关的超像素中心的距离。
%
% It is suggested that use of this function is followed by SPDBSCAN to perform a
% DBSCAN clustering of superpixels.  This results in a simple and fast
% segmentation of an image.建议使用此功能后SPDBSCAN执行超像素的DBSCAN聚类。 这导致图像的简单且快速的分割。
%
% Minor variations from the original algorithm as defined in Achanta et al's
% paper:                     Achanta等人论文中定义的原始算法的微小变化：
%
% - SuperPixel centres are initialised on a hexagonal grid rather than a squareSuperPixel中心在六边形网格而不是方形网格上初始化。 
%   one. This results in a segmentation that will be nominally 6-connected   这导致名义上为6连接的分段，
%   which hopefully facilitates any subsequent post-processing that seeks to 这有希望促进任何后续的寻求合并超像素的后处理。
%   merge superpixels.
% - Initial cluster positions are not shifted to point of lowest gradient初始簇位置不会移动到3x3邻域内
%   within a 3x3 neighbourhood because this will be rendered irrelevant 的最低梯度点，因为这将在第一次
%   first time cluster centres are updated.更新聚类中心时变得无关紧要。
%
% Reference: R. Achanta, A. Shaji, K. Smith, A. Lucchi, P. Fua and
% S. Susstrunk. "SLIC Superpixels Compared to State-of-the-Art Superpixel
% Methods"  PAMI. Vol 34 No 11.  November 2012. pp 2274-2281.
%
% See also: SPDBSCAN, MCLEANUPREGIONS, REGIONADJACENCY, DRAWREGIONBOUNDARIES, RGB2LAB

% Copyright (c) 2013 Peter Kovesi
% www.peterkovesi.com/matlabfns/
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in 
% all copies or substantial portions of the Software.
%
% The Software is provided "as is", without warranty of any kind.

% Feb  2013
% July 2013 Super pixel attributes returned as a structure array

% Note that most of the computation time is not in the clustering, but rather
% in the region cleanup process.请注意，大多数计算时间不在群集中，而是在区域清理过程中。


function [l, Am, Sp, d] = slic(im, k, m, seRadius, colopt, mw, nItr,pa1,pa2)
    
    if ~exist('colopt','var') || isempty(colopt), colopt = 'mean'; end
    if ~exist('mw','var')     || isempty(mw),         mw = 0;      end
    if ~exist('nItr','var')   || isempty(nItr),     nItr = 10;     end
    
    if exist('eim', 'var'), USEDIST = 0; else, USEDIST = 1; end
        
    MEANCENTRE = 1;
    MEDIANCENTRE = 2;
    
    if strcmp(colopt, 'mean')   %strcmp是用于做字符串比较的函数
        centre = MEANCENTRE;
    elseif strcmp(colopt, 'median')
        centre = MEDIANCENTRE;        
    else
        error('Invalid colour centre computation option');%无效的色心计算选项
    end
    
    [rows, cols, chan] = size(im);
    if chan ~= 3
        error('Image must be colour');
    end
    
    % Convert image to L*a*b* colourspace. This gives us a colourspace that is将图像转换为L*a*b*颜色空间。
    % nominally perceptually uniform. This allows us to use the euclidean这为我们提供了一个名义上在感知上 
    % distance between colour coordinates to measure differences between统一的颜色空间。这允许我们使用颜色
    % colours. Note the image becomes double after conversion. We may want to坐标之间的欧氏距离来测量颜色
    % go to signed shorts to save memory.之间的差异。请注意，转换后图像变为双倍。我们可能想签署短以节省内存。
    im = rgb2lab(im); 


%双边滤波%
% pa1 = 15;
% pa2 = 30;
oriL = im(:,:,1);
[width, height]=size(oriL);
sigmaSpatial  = min( width, height ) / pa1;%默认16
samplingSpatial=sigmaSpatial;
sigmaRange = ( max( oriL( : ) ) - min( oriL( : ) ) ) / pa2;%默认10
samplingRange= sigmaRange;
output = bilateralFilter( oriL, oriL, sigmaSpatial, sigmaRange, ...
    samplingSpatial, samplingRange );
oriL = output;

oria = im(:,:,2);
[width, height]=size(oria);
sigmaSpatial  = min( width, height ) / pa1;%默认16
samplingSpatial=sigmaSpatial;
sigmaRange = ( max( oria( : ) ) - min( oria( : ) ) ) / pa2;%默认10
samplingRange= sigmaRange;
output = bilateralFilter( oria, oria, sigmaSpatial, sigmaRange, ...
    samplingSpatial, samplingRange );
oria = output;

orib = im(:,:,3);
[width, height]=size(orib);
sigmaSpatial  = min( width, height ) / pa1;%默认16
samplingSpatial=sigmaSpatial;
sigmaRange = ( max( orib( : ) ) - min( orib( : ) ) ) / pa2;%默认10
samplingRange= sigmaRange;
output = bilateralFilter( orib, orib, sigmaSpatial, sigmaRange, ...
    samplingSpatial, samplingRange );
orib = output;

im=cat(3,oriL,oria,orib);
% im2 = lab2rgb(im);
% figure,imshow(im2,[]);


%     for i=1:rows
%     for j=1:cols
%         rrr = im(i,j,1);
%         r =  double(im(i,j,1));
%         g =  double(im(i,j,2));
%         b =  double(im(i,j,3));
%         sum = double(r+g+b);
%         im(i,j,1) = r/sum*255;
%         im(i,j,2) = g/sum*255;
%         im(i,j,3) = b/sum*255;
%     end
%     end
%     for n = 1:3
%         for i = 1:rows
%             for j = 1:cols
%                 im(i,j,n) = round(im(i,j,n) + 0.1);
%             end
%         end
%     end
%     im = rgb2hsv(im);

%%
tic;
%下采样
[I] = Pyramid(im);
imp = I;
[rowp, colp, chan] = size(imp);
m = m/4;

    % Apply median filtering to colour components if mw has been supplied
    % and/or non-zero如果已提供mw和/或非零，则对颜色分量应用中值滤波
%     if mw
%         if length(mw) == 1
%             mw(2) = mw(1);  % Use same filtering for L and chrominance对L和色度使用相同的过滤
%         end
%         for n = 1:3
%             im(:,:,n) = medfilt2(im(:,:,n), [mw(1) mw(1)]);%中值滤波
%         end
%     end
% tic;
    % Nominal spacing between grid elements assuming hexagonal grid假设六边形网格的网格元素之间的标称间距
    S = sqrt(rowp*colp / (k * sqrt(3)/2));
    
    % Get nodes per row allowing a half column margin at one end that alternates
    % from row to row   每行获取节点允许一端的半列边距从一行到另一行交替
    nodeColp = round(colp/S - 0.5);
    % Given an integer number of nodes per row recompute S  给定每行的整数个节点重新计算S.
    S = colp/(nodeColp + 0.5); 

    % Get number of rows of nodes allowing 0.5 row margin top and bottom获取节点行数，允许0.5行上边距和下边距
    nodeRowp = round(rowp/(sqrt(3)/2*S));
    vSpacing = rowp/nodeRowp;%行采样距离

    % Recompute k 重新计算k
    k = nodeRowp * nodeColp;
    
    
    %正方形采样%
%     S = sqrt(rows*cols / k);
%     nodeCols = round(cols/S);
%     S = cols/(nodeCols);
%     nodeRows = round(rows/S);
%     vSpacing = rows/nodeRows;%行采样距离
%     k = nodeRows * nodeCols;
    
    
    
    % Allocate memory and initialise clusters, labels and distances.分配内存并初始化群集，标签和距离。
    C = zeros(6,k);        % Cluster centre data  1:3 is mean Lab value,  簇中心数据1：3是Lab值，
                           % 4:5 is row, col of centre, 6 is No of pixels 4：5是行，col是中心，6是像素No
    lp = -ones(rowp, colp); % Pixel labels.设标签为-1
    dp = inf(rowp, colp);   % Pixel distances from cluster centres.距离聚类中心的像素距离。
    
    % Initialise clusters on a hexagonal grid在六边形网格上初始化聚类
    kk = 1;%C纵坐标
    r = vSpacing/2;
    
    for ri = 1:nodeRowp
        % Following code alternates the starting column for each row of grid 以下代码交替每行网格点的起始列 
        % points to obtain a hexagonal pattern. Note S and vSpacing are kept 以获得六边形图案。 注意
        % as doubles to prevent errors accumulating across the grid.   S和vSpacing保持双精度，以防止错误累积在网格上。
        if mod(ri,2), c = S/2; else, c = S;  end   %返回数除以除数后的余数，结果与除数具有相同的符号
                                                   %奇数 c = S/2
%         c = S;
        for ci = 1:nodeColp
            cc = round(c); rr = round(r);
            C(1:5, kk) = [squeeze(imp(rr,cc,:)); cc; rr];%squeeze删除单一维度
            c = c+S;
            kk = kk+1;
        end
        
        r = r+vSpacing;
    end
    
    % Now perform the clustering.  10 iterations is suggested but I suspect n
    % could be as small as 2 or even 1   现在执行群集。 建议进行10次迭代，但我怀疑n可能小到2甚至1
    S = round(S);  % We need S to be an integer from now on 从现在开始我们需要S为整数
    

       for kk = 1:k  % for each cluster

           % Get subimage around cluster 获取群集周围的子图像
           rmin = max(C(5,kk)-S, 1);   rmax = min(C(5,kk)+S, rowp); 
           cmin = max(C(4,kk)-S, 1);   cmax = min(C(4,kk)+S, colp); 
           subim = imp(rmin:rmax, cmin:cmax, :);  
           assert(numel(subim) > 0)    %numel返回像素数  
                                   %断言函数assert：在程序中确保某些条件成立，否则调用系统error函数终止运行。
           
           % Compute distances D between C(:,kk) and subimage 计算C（：，kk）和子图像之间的距离D.
           if USEDIST   %判断USEDIST是否为1
               D = dist(C(:, kk), subim, rmin, cmin, S, m);
           else
               D = dist2(C(:, kk), subim, rmin, cmin, S, m, eim, We);
           end

           % If any pixel distance from the cluster centre is less than its 如果距离群集中心的任何像素
           % previous value update its distance and label  距离小于其先前值，则更新其距离和标签
           subd =  dp(rmin:rmax, cmin:cmax);
           subl =  lp(rmin:rmax, cmin:cmax);
           updateMask = D < subd;%先判断大小
           subd(updateMask) = D(updateMask);%subd运算后的最小距离
           subl(updateMask) = kk;
           
           d(rmin:rmax, cmin:cmax) = subd;
           lp(rmin:rmax, cmin:cmax) = subl;%标签           
       end
       
       % Update cluster centres with mean values 使用平均值更新集群中心
       C(:) = 0;%置0
       for r = 1:rowp
           for c = 1:colp
              tmp = [imp(r,c,1); imp(r,c,2); imp(r,c,3); c; r; 1];
              C(:, lp(r,c)) = C(:, lp(r,c)) + tmp;
           end
       end
       
       % Divide by number of pixels in each superpixel to get mean values除以每个超像素中的像素数以获得平均值
       for kk = 1:k 
           C(1:5,kk) = round(C(1:5,kk)/C(6,kk)); 
       end
% toc;
       
       
%%
%正常
m = m*4;
S = S*4;
l = -ones(rows, cols); % Pixel labels.设标签为-1
d = inf(rows, cols);   % Pixel distances from cluster centres.距离聚类中心的像素距离。
% for i = 1:rowp
%     for j = 1:colp
%         if i < rowp && j < colp
%         l(2*i-1,2*j-1) = lp(i,j); 
%         l(2*i,2*j-1) = lp(i,j);
%         l(2*i-1,2*j) = lp(i,j);
%         l(2*i,2*j) = lp(i,j);
%         d(2*i-1,2*j-1) = dp(i,j); 
%         d(2*i,2*j-1) = dp(i,j);
%         d(2*i-1,2*j) = dp(i,j);
%         d(2*i,2*j) = dp(i,j);
%         end
%         if i == rowp && j < colp
%             l(2*i-1,2*j-1) = lp(i,j);
%             l(2*i-1,2*j) = lp(i,j);
%             d(2*i-1,2*j-1) = dp(i,j);
%             d(2*i-1,2*j) = dp(i,j);
%         end
%         if i < rowp && j == colp
%             l(2*i-1,2*j-1) = lp(i,j);
%             l(2*i,2*j-1) = lp(i,j);
%             d(2*i-1,2*j-1) = dp(i,j);
%             d(2*i,2*j-1) = dp(i,j);
%         end
%         if i == rowp && j == colp
%             l(2*i-1,2*j-1) = lp(i,j);
%             d(2*i-1,2*j-1) = dp(i,j);
%         end
%     end
% end
for kk = 1:k
    C(4,kk) = 2*C(4,kk)-1;
    C(5,kk) = 2*C(5,kk)-1;
end

    for n = 1:nItr
       for kk = 1:k  % for each cluster

           % Get subimage around cluster 获取群集周围的子图像
           rmin = max(C(5,kk)-S, 1);   rmax = min(C(5,kk)+S, rows); 
           cmin = max(C(4,kk)-S, 1);   cmax = min(C(4,kk)+S, cols); 
           subim = im(rmin:rmax, cmin:cmax, :);  
           assert(numel(subim) > 0)    %numel返回像素数  
                                   %断言函数assert：在程序中确保某些条件成立，否则调用系统error函数终止运行。
           
           % Compute distances D between C(:,kk) and subimage 计算C（：，kk）和子图像之间的距离D.
           if USEDIST   %判断USEDIST是否为1
               D = dist(C(:, kk), subim, rmin, cmin, S, m);
           else
               D = dist2(C(:, kk), subim, rmin, cmin, S, m, eim, We);
           end

           % If any pixel distance from the cluster centre is less than its 如果距离群集中心的任何像素
           % previous value update its distance and label  距离小于其先前值，则更新其距离和标签
           subd =  d(rmin:rmax, cmin:cmax);
           subl =  l(rmin:rmax, cmin:cmax);
           updateMask = D < subd;%先判断大小
           subd(updateMask) = D(updateMask);%subd运算后的最小距离
           subl(updateMask) = kk;
           
           d(rmin:rmax, cmin:cmax) = subd;
           l(rmin:rmax, cmin:cmax) = subl;%标签           
       end
       
       % Update cluster centres with mean values 使用平均值更新集群中心
       C(:) = 0;%置0
       for r = 1:rows
           for c = 1:cols
              tmp = [im(r,c,1); im(r,c,2); im(r,c,3); c; r; 1];
              C(:, l(r,c)) = C(:, l(r,c)) + tmp;
           end
       end
       
       % Divide by number of pixels in each superpixel to get mean values除以每个超像素中的像素数以获得平均值
       for kk = 1:k 
           C(1:5,kk) = round(C(1:5,kk)/C(6,kk)); 
       end
       
       % Note the residual error, E, is not calculated because we are using a
       % fixed number of iterations 请注意，不计算残差E，因为我们使用的是固定次数的迭代
    end
 %%      

    
    % Cleanup small orphaned regions and 'spurs' on each region using 使用每个标记区域上的形态开口清理每个 
    % morphological opening on each labeled region.  The cleaned up regions are 区域上的小孤立区域和“刺”。 
    % assigned to the nearest cluster. The regions are renumbered and the  已清理的区域将分配给最近的群集。 
    % adjacency matrix regenerated.  This is needed because the cleanup is 区域重新编号并重新生成邻接矩阵。
    % likely to change the number of labeled regions.这是必需的，因为清理可能会改变标记区域的数量。
    if seRadius %1执行
        [l, Am] = mcleanupregions(l, seRadius);%l标签
    else
        l = makeregionsdistinct(l);
        [l, minLabel, maxLabel] = renumberregions(l);
        Am = regionadjacency(l);    
    end

    % Recompute the final superpixel attributes and write information into
    % the Sp struct array.  重新计算最终的超像素属性并将信息写入Sp结构数组。
    N = length(Am);
    Sp = struct('L', cell(1,N), 'a', cell(1,N), 'b', cell(1,N), ...
                'stdL', cell(1,N), 'stda', cell(1,N), 'stdb', cell(1,N), ...
                'r', cell(1,N), 'c', cell(1,N), 'N', cell(1,N));
    [X,Y] = meshgrid(1:cols, 1:rows);%用于生成网格采样点的函数
    L = im(:,:,1);    
    A = im(:,:,2);    
    B = im(:,:,3);    
    for n = 1:N
        mask = l==n;%先运行==
        nm = sum(mask(:));
        if centre == MEANCENTRE     
            Sp(n).L = sum(L(mask))/nm;
            Sp(n).a = sum(A(mask))/nm;
            Sp(n).b = sum(B(mask))/nm;
            
        elseif centre == MEDIANCENTRE
            Sp(n).L = median(L(mask));
            Sp(n).a = median(A(mask));
            Sp(n).b = median(B(mask));
        end
        
        Sp(n).r = sum(Y(mask))/nm;
        Sp(n).c = sum(X(mask))/nm;
        
        % Compute standard deviations of the colour components of each super计算每个超像素的颜色分量的标准偏差。  
        % pixel. This can be used by code seeking to merge superpixels into这可以由寻求将超像素合并到图像片段中的代码使用。
        % image segments.  Note these are calculated relative to the mean colour注意，这些是相对于平均颜色分量计算的，
        % component irrespective of the centre being calculated from the mean or而与从平均或中值颜色分量值计算的中心无关。
        % median colour component values.
        Sp(n).stdL = std(L(mask));
        Sp(n).stda = std(A(mask));
        Sp(n).stdb = std(B(mask));

        Sp(n).N = nm;  % Record number of pixels in superpixel too.记录超像素中的像素数。
    end
    
% %-- dist -------------------------------------------
% %
% % Usage:  D = dist(C, im, r1, c1, S, m)
% % 
% % Arguments:   C - Cluster being considered正在考虑集群
% %             im - sub-image surrounding cluster centre围绕集群中心的子图像
% %         r1, c1 - row and column of top left corner of sub image within the
% %                  overall image.   整个图像中子图像左上角的行和列。
% %              S - grid spacing 网格间距
% %              m - weighting factor between colour and spatial differences.颜色和空间差异之间的加权因子。
% %
% % Returns:     D - Distance image giving distance of every pixel in the
% %                  subimage from the cluster centre  距离图像给出子图像中每个像素距离聚类中心的距离
% %
% % Distance = sqrt( dc^2 + (ds/S)^2*m^2 )
% % where:
% % dc = sqrt(dl^2 + da^2 + db^2)  % Colour distance
% % ds = sqrt(dx^2 + dy^2)         % Spatial distance
% %
% % m is a weighting factor representing the nominal maximum colour distancem是表示预期的标称最大颜色距离
% % expected so that one can rank colour similarity relative to distance的加权因子，因此可以相对于距离相似性对颜色相似性
% % similarity.  try m in the range [1-40] for L*a*b* space进行排序。 对于L * a * b *空间，尝试在[1-40]范围内的m
% %
% % ?? Might be worth trying the Geometric Mean instead ?? 可能值得尝试几何平均值
% %  Distance = sqrt(dc * ds)
% % but having a factor 'm' to play with is probably handy  但有一个因素'm'可能很方便
% 
% % This code could be more efficient  此代码可能更有效
% 
% function D = dist(C, im, r1, c1, S, m)
% 
%     % Squared spatial distance  平方的空间距离
%     %    ds is a fixed 'image' we should be able to exploit this  ds是一个固定的'图像'，我们应该能够
%     %    and use a fixed meshgrid for much of the time somehow...  利用它并在某种程度上使用固定的网格...
%     [rows, cols, chan] = size(im);
%     [x,y] = meshgrid(c1:(c1+cols-1), r1:(r1+rows-1));
%     x = x-C(4);  % x and y dist from cluster centre
%     y = y-C(5);
%     ds2 = x.^2 + y.^2;
%     
%     % Squared colour difference平方色差
%     for n = 1:3
%         im(:,:,n) = (im(:,:,n)-C(n)).^2;
%     end
%     dc2 = sum(im,3);
%     
%     D = sqrt(dc2 + ds2/S^2*m^2);
%     
%     
%     
% %--- dist2 ------------------------------------------
% %
% % Usage:  D = dist2(C, im, r1, c1, S, m, eim)
% % 
% % Arguments:   C - Cluster being considered
% %             im - sub-image surrounding cluster centre
% %         r1, c1 - row and column of top left corner of sub image within the
% %                  overall image.
% %              S - grid spacing
% %              m - weighting factor between colour and spatial differences.
% %            eim - Edge strength sub-image corresponding to im 边缘强度子图像对应于im
% %
% % Returns:     D - Distance image giving distance of every pixel in the
% %                  subimage from the cluster centre
% %
% % Distance = sqrt( dc^2 + (ds/S)^2*m^2 )
% % where:
% % dc = sqrt(dl^2 + da^2 + db^2)  % Colour distance
% % ds = sqrt(dx^2 + dy^2)         % Spatial distance
% %
% % m is a weighting factor representing the nominal maximum colour distance
% % expected so that one can rank colour similarity relative to distance
% % similarity.  try m in the range [1-40] for L*a*b* space
% %
% 
% function D = dist2(C, im, r1, c1, S, m, eim, We)
% 
%     % Squared spatial distance
%     %    ds is a fixed 'image' we should be able to exploit this
%     %    and use a fixed meshgrid for much of the time somehow...
%     [rows, cols, chan] = size(im);
%     [x,y] = meshgrid(c1:(c1+cols-1), r1:(r1+rows-1));
%     x = x-C(4);
%     y = y-C(5);
%     ds2 = x.^2 + y.^2;
%     
%     % Squared colour difference
%     for n = 1:3
%         im(:,:,n) = (im(:,:,n)-C(n)).^2;
%     end
%     dc2 = sum(im,3);
%     
%     % Combine colour and spatial distance measure结合颜色和空间距离测量
%     D = sqrt(dc2 + ds2/S^2*m^2);
%     
%     % for every pixel in the subimage call improfile to the cluster centre对于子图像调用中的每个像素，
%     % and use the largest value as the 'edge distance' 对集群中心进行格式化，并使用最大值作为“边距”
%     rCentre = C(5)-r1;   % Cluster centre coords relative to this sub-image 聚类中心相对于该子图像的坐标
%     cCentre = C(4)-c1;
%     de = zeros(rows,cols);
%     for r = 1:rows
%         for c = 1:cols
%             v = improfile(eim,[c cCentre], [r rCentre]);
%             de(r,c) = max(v);
%         end
%     end
% 
%     % Combine edge distance with weight, We with total Distance.将边距与重量相结合，我们与总距离相结合。
%     D = D + We * de;
    
