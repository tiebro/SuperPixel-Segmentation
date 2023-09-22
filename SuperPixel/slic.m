% SLIC Simple Linear Iterative Clustering SuperPixels�����Ե������೬����
%
% Implementation of Achanta, Shaji, Smith, Lucchi, Fua and Susstrunk's
% SLIC Superpixels
%
% Usage:   [l, Am, Sp, d] = slic(im, k, m, seRadius, colopt, mw)
%
% Arguments:  im - Image to be segmented.Ҫ�ָ��ͼ��
%              k - Number of desired superpixels. Note that this is nominal���賬���ص�������ע�⣬���Ǳ�Ƶģ�
%                  the actual number of superpixels generated will generally�����ɵĳ����ص�ʵ������ͨ�����Դ�
%                  be a bit larger, especially if parameter m is small.�ر����������m��С��
%              m - Weighting factor between colour and spatial   ��ɫ�Ϳռ����֮��ļ�Ȩ���ӡ�   
%                  differences. Values from about 5 to 40 are useful.  Use a  ��Լ5��40��ֵ�����õġ�ʹ�ý�
%                  large value to enforce superpixels with more regular and���ֵ��ǿ��ʹ�ø�����͸�ƽ����
%                  smoother shapes. Try a value of 10 to start with.��״�ĳ����ء�����ʹ��ֵ10��ʼ��
%       seRadius - Regions morphologically smaller than this are merged with  ��̬ѧ��С�ڴ˵���������������ϲ��� 
%                  adjacent regions. Try a value of 1 or 1.5.  Use 0 to       ����ֵ1��1.5�� ʹ��0���á�
%                  disable.
%         colopt - String 'mean' or 'median' indicating how the cluster  �ַ���'mean'��'median'��ʾӦ���
%                  colour centre should be computed. Defaults to 'mean'  �������ɫ���ġ� Ĭ��Ϊ'��ֵ'
%             mw - Optional median filtering window size.  Image compression ��ѡ����ֵ���˴��ڴ�С��ͼ��ѹ��    
%                  can result in noticeable artifacts in the a*b* components ���ܻᵼ��ͼ���a* b*����������Ե�
%                  of the image.  Median filtering can reduce this. mw can be 覴á���ֵ���˿��Լ������������
%                  a single value in which case the same median filtering ismw �����ǵ���ֵ������������£�
%                  applied to each L* a* and b* components.  Alternatively it ��ÿ��L*a*��b*����Ӧ����ͬ����ֵ
%                  can be a 2-vector where mw(1) specifies the median  �˲������ߣ���������2����������mw��1��
%                  filtering window to be applied to L* and mw(2) is the ָ��ҪӦ����L*����ֵ�˲����ڣ�
%                  median filtering window to be applied to a* and b*. ��mw��2����ҪӦ����a*��b*����ֵ�˲����ڡ�
%
% Returns:     l - Labeled image of superpixels. Labels range from 1 to k.��ǩͼ��ĳ����ء� ��ǩ��Χ��1��k��
%             Am - Adjacency matrix of segments.  Am(i, j) indicates whether  �ָ���ڽӾ��� Am��i��j��
%                  segments labeled i and j are connected/adjacent      ��ʾ���Ϊi��j�Ķ��Ƿ�����/����
%             Sp - Superpixel attribute structure array with fields:  �����ֶεĳ��������Խṹ���飺
%                   L  - Mean L* value ƽ��L*ֵ
%                   a  - Mean a* value ƽ��a*ֵ
%                   b  - Mean b* value ƽ��b*ֵ
%                   r  - Mean row value ƽ��rowֵ
%                   c  - Mean column value ƽ��columnֵ
%                   stdL  - Standard deviation of L*     L*�ı�׼ƫ��
%                   stda  - Standard deviation of a*     a*�ı�׼ƫ�� 
%                   stdb  - Standard deviation of b*     b*�ı�׼ƫ��
%                   N - Number of pixels
%                   edges - List of edge numbers that bound each  ��ÿ�������صı����б� 
%                           superpixel. This field is allocated, but not set, ���ֶ���SLIC���䣬
%                           by SLIC. Use  for this.  ��δ���á�ʹ��SPEDGES��
%              d - Distance image giving the distance each pixel is from its
%                  associated superpixel centre.����ͼ�����ÿ�����ؾ�����صĳ��������ĵľ��롣
%
% It is suggested that use of this function is followed by SPDBSCAN to perform a
% DBSCAN clustering of superpixels.  This results in a simple and fast
% segmentation of an image.����ʹ�ô˹��ܺ�SPDBSCANִ�г����ص�DBSCAN���ࡣ �⵼��ͼ��ļ��ҿ��ٵķָ
%
% Minor variations from the original algorithm as defined in Achanta et al's
% paper:                     Achanta���������ж����ԭʼ�㷨��΢С�仯��
%
% - SuperPixel centres are initialised on a hexagonal grid rather than a squareSuperPixel��������������������Ƿ��������ϳ�ʼ���� 
%   one. This results in a segmentation that will be nominally 6-connected   �⵼��������Ϊ6���ӵķֶΣ�
%   which hopefully facilitates any subsequent post-processing that seeks to ����ϣ���ٽ��κκ�����Ѱ��ϲ������صĺ���
%   merge superpixels.
% - Initial cluster positions are not shifted to point of lowest gradient��ʼ��λ�ò����ƶ���3x3������
%   within a 3x3 neighbourhood because this will be rendered irrelevant ������ݶȵ㣬��Ϊ�⽫�ڵ�һ��
%   first time cluster centres are updated.���¾�������ʱ����޹ؽ�Ҫ��
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
% in the region cleanup process.��ע�⣬���������ʱ�䲻��Ⱥ���У�������������������С�


function [l, Am, Sp, d] = slic(im, k, m, seRadius, colopt, mw, nItr,pa1,pa2)
    
    if ~exist('colopt','var') || isempty(colopt), colopt = 'mean'; end
    if ~exist('mw','var')     || isempty(mw),         mw = 0;      end
    if ~exist('nItr','var')   || isempty(nItr),     nItr = 10;     end
    
    if exist('eim', 'var'), USEDIST = 0; else, USEDIST = 1; end
        
    MEANCENTRE = 1;
    MEDIANCENTRE = 2;
    
    if strcmp(colopt, 'mean')   %strcmp���������ַ����Ƚϵĺ���
        centre = MEANCENTRE;
    elseif strcmp(colopt, 'median')
        centre = MEDIANCENTRE;        
    else
        error('Invalid colour centre computation option');%��Ч��ɫ�ļ���ѡ��
    end
    
    [rows, cols, chan] = size(im);
    if chan ~= 3
        error('Image must be colour');
    end
    
    % Convert image to L*a*b* colourspace. This gives us a colourspace that is��ͼ��ת��ΪL*a*b*��ɫ�ռ䡣
    % nominally perceptually uniform. This allows us to use the euclidean��Ϊ�����ṩ��һ���������ڸ�֪�� 
    % distance between colour coordinates to measure differences betweenͳһ����ɫ�ռ䡣����������ʹ����ɫ
    % colours. Note the image becomes double after conversion. We may want to����֮���ŷ�Ͼ�����������ɫ
    % go to signed shorts to save memory.֮��Ĳ��졣��ע�⣬ת����ͼ���Ϊ˫�������ǿ�����ǩ����Խ�ʡ�ڴ档
    im = rgb2lab(im); 


%˫���˲�%
% pa1 = 15;
% pa2 = 30;
oriL = im(:,:,1);
[width, height]=size(oriL);
sigmaSpatial  = min( width, height ) / pa1;%Ĭ��16
samplingSpatial=sigmaSpatial;
sigmaRange = ( max( oriL( : ) ) - min( oriL( : ) ) ) / pa2;%Ĭ��10
samplingRange= sigmaRange;
output = bilateralFilter( oriL, oriL, sigmaSpatial, sigmaRange, ...
    samplingSpatial, samplingRange );
oriL = output;

oria = im(:,:,2);
[width, height]=size(oria);
sigmaSpatial  = min( width, height ) / pa1;%Ĭ��16
samplingSpatial=sigmaSpatial;
sigmaRange = ( max( oria( : ) ) - min( oria( : ) ) ) / pa2;%Ĭ��10
samplingRange= sigmaRange;
output = bilateralFilter( oria, oria, sigmaSpatial, sigmaRange, ...
    samplingSpatial, samplingRange );
oria = output;

orib = im(:,:,3);
[width, height]=size(orib);
sigmaSpatial  = min( width, height ) / pa1;%Ĭ��16
samplingSpatial=sigmaSpatial;
sigmaRange = ( max( orib( : ) ) - min( orib( : ) ) ) / pa2;%Ĭ��10
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
%�²���
[I] = Pyramid(im);
imp = I;
[rowp, colp, chan] = size(imp);
m = m/4;

    % Apply median filtering to colour components if mw has been supplied
    % and/or non-zero������ṩmw��/����㣬�����ɫ����Ӧ����ֵ�˲�
%     if mw
%         if length(mw) == 1
%             mw(2) = mw(1);  % Use same filtering for L and chrominance��L��ɫ��ʹ����ͬ�Ĺ���
%         end
%         for n = 1:3
%             im(:,:,n) = medfilt2(im(:,:,n), [mw(1) mw(1)]);%��ֵ�˲�
%         end
%     end
% tic;
    % Nominal spacing between grid elements assuming hexagonal grid�������������������Ԫ��֮��ı�Ƽ��
    S = sqrt(rowp*colp / (k * sqrt(3)/2));
    
    % Get nodes per row allowing a half column margin at one end that alternates
    % from row to row   ÿ�л�ȡ�ڵ�����һ�˵İ��б߾��һ�е���һ�н���
    nodeColp = round(colp/S - 0.5);
    % Given an integer number of nodes per row recompute S  ����ÿ�е��������ڵ����¼���S.
    S = colp/(nodeColp + 0.5); 

    % Get number of rows of nodes allowing 0.5 row margin top and bottom��ȡ�ڵ�����������0.5���ϱ߾���±߾�
    nodeRowp = round(rowp/(sqrt(3)/2*S));
    vSpacing = rowp/nodeRowp;%�в�������

    % Recompute k ���¼���k
    k = nodeRowp * nodeColp;
    
    
    %�����β���%
%     S = sqrt(rows*cols / k);
%     nodeCols = round(cols/S);
%     S = cols/(nodeCols);
%     nodeRows = round(rows/S);
%     vSpacing = rows/nodeRows;%�в�������
%     k = nodeRows * nodeCols;
    
    
    
    % Allocate memory and initialise clusters, labels and distances.�����ڴ沢��ʼ��Ⱥ������ǩ�;��롣
    C = zeros(6,k);        % Cluster centre data  1:3 is mean Lab value,  ����������1��3��Labֵ��
                           % 4:5 is row, col of centre, 6 is No of pixels 4��5���У�col�����ģ�6������No
    lp = -ones(rowp, colp); % Pixel labels.���ǩΪ-1
    dp = inf(rowp, colp);   % Pixel distances from cluster centres.����������ĵ����ؾ��롣
    
    % Initialise clusters on a hexagonal grid�������������ϳ�ʼ������
    kk = 1;%C������
    r = vSpacing/2;
    
    for ri = 1:nodeRowp
        % Following code alternates the starting column for each row of grid ���´��뽻��ÿ����������ʼ�� 
        % points to obtain a hexagonal pattern. Note S and vSpacing are kept �Ի��������ͼ���� ע��
        % as doubles to prevent errors accumulating across the grid.   S��vSpacing����˫���ȣ��Է�ֹ�����ۻ��������ϡ�
        if mod(ri,2), c = S/2; else, c = S;  end   %���������Գ��������������������������ͬ�ķ���
                                                   %���� c = S/2
%         c = S;
        for ci = 1:nodeColp
            cc = round(c); rr = round(r);
            C(1:5, kk) = [squeeze(imp(rr,cc,:)); cc; rr];%squeezeɾ����һά��
            c = c+S;
            kk = kk+1;
        end
        
        r = r+vSpacing;
    end
    
    % Now perform the clustering.  10 iterations is suggested but I suspect n
    % could be as small as 2 or even 1   ����ִ��Ⱥ���� �������10�ε��������һ���n����С��2����1
    S = round(S);  % We need S to be an integer from now on �����ڿ�ʼ������ҪSΪ����
    

       for kk = 1:k  % for each cluster

           % Get subimage around cluster ��ȡȺ����Χ����ͼ��
           rmin = max(C(5,kk)-S, 1);   rmax = min(C(5,kk)+S, rowp); 
           cmin = max(C(4,kk)-S, 1);   cmax = min(C(4,kk)+S, colp); 
           subim = imp(rmin:rmax, cmin:cmax, :);  
           assert(numel(subim) > 0)    %numel����������  
                                   %���Ժ���assert���ڳ�����ȷ��ĳЩ�����������������ϵͳerror������ֹ���С�
           
           % Compute distances D between C(:,kk) and subimage ����C������kk������ͼ��֮��ľ���D.
           if USEDIST   %�ж�USEDIST�Ƿ�Ϊ1
               D = dist(C(:, kk), subim, rmin, cmin, S, m);
           else
               D = dist2(C(:, kk), subim, rmin, cmin, S, m, eim, We);
           end

           % If any pixel distance from the cluster centre is less than its �������Ⱥ�����ĵ��κ�����
           % previous value update its distance and label  ����С������ǰֵ������������ͱ�ǩ
           subd =  dp(rmin:rmax, cmin:cmax);
           subl =  lp(rmin:rmax, cmin:cmax);
           updateMask = D < subd;%���жϴ�С
           subd(updateMask) = D(updateMask);%subd��������С����
           subl(updateMask) = kk;
           
           d(rmin:rmax, cmin:cmax) = subd;
           lp(rmin:rmax, cmin:cmax) = subl;%��ǩ           
       end
       
       % Update cluster centres with mean values ʹ��ƽ��ֵ���¼�Ⱥ����
       C(:) = 0;%��0
       for r = 1:rowp
           for c = 1:colp
              tmp = [imp(r,c,1); imp(r,c,2); imp(r,c,3); c; r; 1];
              C(:, lp(r,c)) = C(:, lp(r,c)) + tmp;
           end
       end
       
       % Divide by number of pixels in each superpixel to get mean values����ÿ���������е��������Ի��ƽ��ֵ
       for kk = 1:k 
           C(1:5,kk) = round(C(1:5,kk)/C(6,kk)); 
       end
% toc;
       
       
%%
%����
m = m*4;
S = S*4;
l = -ones(rows, cols); % Pixel labels.���ǩΪ-1
d = inf(rows, cols);   % Pixel distances from cluster centres.����������ĵ����ؾ��롣
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

           % Get subimage around cluster ��ȡȺ����Χ����ͼ��
           rmin = max(C(5,kk)-S, 1);   rmax = min(C(5,kk)+S, rows); 
           cmin = max(C(4,kk)-S, 1);   cmax = min(C(4,kk)+S, cols); 
           subim = im(rmin:rmax, cmin:cmax, :);  
           assert(numel(subim) > 0)    %numel����������  
                                   %���Ժ���assert���ڳ�����ȷ��ĳЩ�����������������ϵͳerror������ֹ���С�
           
           % Compute distances D between C(:,kk) and subimage ����C������kk������ͼ��֮��ľ���D.
           if USEDIST   %�ж�USEDIST�Ƿ�Ϊ1
               D = dist(C(:, kk), subim, rmin, cmin, S, m);
           else
               D = dist2(C(:, kk), subim, rmin, cmin, S, m, eim, We);
           end

           % If any pixel distance from the cluster centre is less than its �������Ⱥ�����ĵ��κ�����
           % previous value update its distance and label  ����С������ǰֵ������������ͱ�ǩ
           subd =  d(rmin:rmax, cmin:cmax);
           subl =  l(rmin:rmax, cmin:cmax);
           updateMask = D < subd;%���жϴ�С
           subd(updateMask) = D(updateMask);%subd��������С����
           subl(updateMask) = kk;
           
           d(rmin:rmax, cmin:cmax) = subd;
           l(rmin:rmax, cmin:cmax) = subl;%��ǩ           
       end
       
       % Update cluster centres with mean values ʹ��ƽ��ֵ���¼�Ⱥ����
       C(:) = 0;%��0
       for r = 1:rows
           for c = 1:cols
              tmp = [im(r,c,1); im(r,c,2); im(r,c,3); c; r; 1];
              C(:, l(r,c)) = C(:, l(r,c)) + tmp;
           end
       end
       
       % Divide by number of pixels in each superpixel to get mean values����ÿ���������е��������Ի��ƽ��ֵ
       for kk = 1:k 
           C(1:5,kk) = round(C(1:5,kk)/C(6,kk)); 
       end
       
       % Note the residual error, E, is not calculated because we are using a
       % fixed number of iterations ��ע�⣬������в�E����Ϊ����ʹ�õ��ǹ̶������ĵ���
    end
 %%      

    
    % Cleanup small orphaned regions and 'spurs' on each region using ʹ��ÿ����������ϵ���̬��������ÿ�� 
    % morphological opening on each labeled region.  The cleaned up regions are �����ϵ�С��������͡��̡��� 
    % assigned to the nearest cluster. The regions are renumbered and the  ����������򽫷���������Ⱥ���� 
    % adjacency matrix regenerated.  This is needed because the cleanup is �������±�Ų����������ڽӾ���
    % likely to change the number of labeled regions.���Ǳ���ģ���Ϊ������ܻ�ı��������������
    if seRadius %1ִ��
        [l, Am] = mcleanupregions(l, seRadius);%l��ǩ
    else
        l = makeregionsdistinct(l);
        [l, minLabel, maxLabel] = renumberregions(l);
        Am = regionadjacency(l);    
    end

    % Recompute the final superpixel attributes and write information into
    % the Sp struct array.  ���¼������յĳ��������Բ�����Ϣд��Sp�ṹ���顣
    N = length(Am);
    Sp = struct('L', cell(1,N), 'a', cell(1,N), 'b', cell(1,N), ...
                'stdL', cell(1,N), 'stda', cell(1,N), 'stdb', cell(1,N), ...
                'r', cell(1,N), 'c', cell(1,N), 'N', cell(1,N));
    [X,Y] = meshgrid(1:cols, 1:rows);%�����������������ĺ���
    L = im(:,:,1);    
    A = im(:,:,2);    
    B = im(:,:,3);    
    for n = 1:N
        mask = l==n;%������==
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
        
        % Compute standard deviations of the colour components of each super����ÿ�������ص���ɫ�����ı�׼ƫ�  
        % pixel. This can be used by code seeking to merge superpixels into�������Ѱ�󽫳����غϲ���ͼ��Ƭ���еĴ���ʹ�á�
        % image segments.  Note these are calculated relative to the mean colourע�⣬��Щ�������ƽ����ɫ��������ģ�
        % component irrespective of the centre being calculated from the mean or�����ƽ������ֵ��ɫ����ֵ����������޹ء�
        % median colour component values.
        Sp(n).stdL = std(L(mask));
        Sp(n).stda = std(A(mask));
        Sp(n).stdb = std(B(mask));

        Sp(n).N = nm;  % Record number of pixels in superpixel too.��¼�������е���������
    end
    
% %-- dist -------------------------------------------
% %
% % Usage:  D = dist(C, im, r1, c1, S, m)
% % 
% % Arguments:   C - Cluster being considered���ڿ��Ǽ�Ⱥ
% %             im - sub-image surrounding cluster centreΧ�Ƽ�Ⱥ���ĵ���ͼ��
% %         r1, c1 - row and column of top left corner of sub image within the
% %                  overall image.   ����ͼ������ͼ�����Ͻǵ��к��С�
% %              S - grid spacing ������
% %              m - weighting factor between colour and spatial differences.��ɫ�Ϳռ����֮��ļ�Ȩ���ӡ�
% %
% % Returns:     D - Distance image giving distance of every pixel in the
% %                  subimage from the cluster centre  ����ͼ�������ͼ����ÿ�����ؾ���������ĵľ���
% %
% % Distance = sqrt( dc^2 + (ds/S)^2*m^2 )
% % where:
% % dc = sqrt(dl^2 + da^2 + db^2)  % Colour distance
% % ds = sqrt(dx^2 + dy^2)         % Spatial distance
% %
% % m is a weighting factor representing the nominal maximum colour distancem�Ǳ�ʾԤ�ڵı�������ɫ����
% % expected so that one can rank colour similarity relative to distance�ļ�Ȩ���ӣ���˿�������ھ��������Զ���ɫ������
% % similarity.  try m in the range [1-40] for L*a*b* space�������� ����L * a * b *�ռ䣬������[1-40]��Χ�ڵ�m
% %
% % ?? Might be worth trying the Geometric Mean instead ?? ����ֵ�ó��Լ���ƽ��ֵ
% %  Distance = sqrt(dc * ds)
% % but having a factor 'm' to play with is probably handy  ����һ������'m'���ܷܺ���
% 
% % This code could be more efficient  �˴�����ܸ���Ч
% 
% function D = dist(C, im, r1, c1, S, m)
% 
%     % Squared spatial distance  ƽ���Ŀռ����
%     %    ds is a fixed 'image' we should be able to exploit this  ds��һ���̶���'ͼ��'������Ӧ���ܹ�
%     %    and use a fixed meshgrid for much of the time somehow...  ����������ĳ�̶ֳ���ʹ�ù̶�������...
%     [rows, cols, chan] = size(im);
%     [x,y] = meshgrid(c1:(c1+cols-1), r1:(r1+rows-1));
%     x = x-C(4);  % x and y dist from cluster centre
%     y = y-C(5);
%     ds2 = x.^2 + y.^2;
%     
%     % Squared colour differenceƽ��ɫ��
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
% %            eim - Edge strength sub-image corresponding to im ��Եǿ����ͼ���Ӧ��im
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
%     % Combine colour and spatial distance measure�����ɫ�Ϳռ�������
%     D = sqrt(dc2 + ds2/S^2*m^2);
%     
%     % for every pixel in the subimage call improfile to the cluster centre������ͼ������е�ÿ�����أ�
%     % and use the largest value as the 'edge distance' �Լ�Ⱥ���Ľ��и�ʽ������ʹ�����ֵ��Ϊ���߾ࡱ
%     rCentre = C(5)-r1;   % Cluster centre coords relative to this sub-image ������������ڸ���ͼ�������
%     cCentre = C(4)-c1;
%     de = zeros(rows,cols);
%     for r = 1:rows
%         for c = 1:cols
%             v = improfile(eim,[c cCentre], [r rCentre]);
%             de(r,c) = max(v);
%         end
%     end
% 
%     % Combine edge distance with weight, We with total Distance.���߾����������ϣ��������ܾ������ϡ�
%     D = D + We * de;
    
