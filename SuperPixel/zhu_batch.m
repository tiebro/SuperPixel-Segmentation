clear;clc;close all;
% [fileName,pathName] = uigetfile('*.*','Please select an image');%文件筐，选择文件
% if(fileName)
%     fileName = strcat(pathName,fileName);
%     fileName = lower(fileName);%一致的小写字母形式
% else 
%     msgbox('Please select an image');
%     return; %退出程序
% end
% im = imread(fileName);
% tic;

img1 = 'C:\Users\15650\Desktop\颈椎相关资料\超像素算法分割颈椎曲率计算\MRI\MRI_data\上海市中医院\上海市中医院-双边滤波\';
% img1 = 'C:\Users\15650\Desktop\颈椎相关资料\超像素算法分割颈椎曲率计算\MRI\MRI_data\上海市中医院\1\';
img_Dir = dir(strcat(img1,'*.bmp'));
% mask = 'C:\Users\15650\Desktop\颈椎相关资料\超像素算法分割颈椎曲率计算\MRI\超像素分割 - 4 双边+金字塔+六边形\mask\';
% mask_Dir = dir(strcat(mask,'*.jpg'));
t = 1;
for pn = 1:length(img_Dir)
    name1 = img_Dir(pn).name;
    im = imread(strcat(img1,name1));
    fprintf('%d %s\n',pn,strcat(img1,name1));% 显示正在处理的图像名
    im_name = char(name1);
    input_age = string(im_name(:,8:9));
    input_sex = string(im_name(:,10));
%     name1='裁剪23-T1-6 3';%裁剪2 图片5  裁剪23-T1-6 3  裁剪1135-T1-7 3
%     name2 = '真实裁剪23-T1-6 3-2';
    % tic;
%     im=imread([name1,'.jpg']);%jpg png
    % figure,imshow(im);
    gray = rgb2gray(im);

    k = 200;%初始超像素总数
    m = 10;%颜色和空间差异之间的加权因子，5-20
    seRadius = 1;%小于此值的区域与相邻区域合并，seRadius = 1.5，可以设置为0
    colopt = 'mean';%计算簇中心分方式，可用'mean'或'median'
    mw = 3;%均值/中值过滤窗口大小
    nItr = 3;%迭代次数
    pa1 = 9;%pa2足够小时，越大越平滑
    pa2 = 30;
    eim = 0;%对距离公式有影响,暂时用不上
    op1 = 1;%op1为曲度测量方法1
    op2 = 1;%op2为曲度测量方法2
    Ec = 9; %5-15



    %直方图均衡化图像增强
    % im_R=im(:,:,1);
    % im_G= im(:,:,2);
    % im_B= im(:,:,3);
    % r=histeq(im_R);%对个分量直方图均衡化，得到个分量均衡化图像?
    % g=histeq(im_G);
    % b=histeq(im_B);
    % im = cat(3,r,g,b);
    % figure,imshow(im);
    % figure,imshow(r);
    % figure,imshow(g);
    % figure,imshow(b);
    % figure,
    % subplot(131),imhist(r,64);
    % subplot(132),imhist(g,64);
    % subplot(133),imhist(b,64);


    % im = im2double(im); % 读入图像
    % h = [1,1,1; 1,-8,1; 1,1,1]; % 8拉普拉斯算子
    % % h = [0,1,0; 1,-4,1; 0,1,0]; % 4拉普拉斯算子
    % im_R=im(:,:,1);
    % im_G= im(:,:,2);
    % im_B= im(:,:,3);
    % J_R = conv2(im_R, h, 'same');  % 通过卷积对图像进行锐化滤波
    % J_G = conv2(im_G, h, 'same');  % 通过卷积对图像进行锐化滤波
    % J_B = conv2(im_B, h, 'same');  % 通过卷积对图像进行锐化滤波
    % J = cat(3,J_R,J_G,J_B);
    % im = im - J;
    % figure,imshow(im,[]);


    [l, Am, Sp, d] = slic(im, k, m, seRadius, colopt, mw, nItr,pa1,pa2);%标签图像的超像素。标签范围从1到k。
    % figure,imshow(drawregionboundaries(l,im,[255 255 255]));
    % DisplaySuperpixel(l,im,name1);%轮廓较粗
    % super = drawregionboundaries(l,im,[255 255 255]);
    % imwrite(super,'G:\MRI\图片\our\k = 700 自动分割 金\自动分割1221-T1-7 3-1.jpg');
    toc;

    E3 = 2.1;
    [lc, C, regionsC] = spdbscan(l, Sp, Am, Ec, gray, E3);%新标记的图像对应于新的超像素聚集区域
    figure,imshow(drawregionboundaries(lc,im,[255 255 255]));
    % hebing = drawregionboundaries(lc,im,[255 255 255]);
    % imwrite(hebing,'G:\MRI\论文材料2\实验二\our\k = 400 自动分割 金\合并23-T1-6 3-1.jpg');




    % regionadjacency 4联通
    pos1 = [];%放分割中心点坐标
    [x,y] = ginput;
    [row,col] = size(im(:,:,1));
    imzhuiti = [];
    imzhuiti = im;
    img = zeros(row,col);
    for n = 1:3
        img(:,:,n) = zeros(row,col);
    end
    a = length(x);
    for i = 1:a
        X = round(x(i,1));Y = round(y(i,1));
        b = lc(Y,X);
        [u,v]=find(lc==b);
        d = length(u);
        pos1(i,1) = round(mean(u));%纵坐标
        pos1(i,2) = round(mean(v));%横坐标
        for j = 1:d
            for n = 1:3
            img(u(j,1),v(j,1),n) = im(u(j,1),v(j,1),n);
            end
        end
        for j = 1:d                                                                                                                                                         
            imzhuiti(u(j,1),v(j,1),1) = 255;
            imzhuiti(u(j,1),v(j,1),2) = 0;
            imzhuiti(u(j,1),v(j,1),3) = 0;
        end
    end

    % figure,imshow(img,[]);%分割的黑白图像
    %为了第四章 k=200
    % imwrite(img,'H:\MRI\论文材料2\第四章实验\自动分割1135-T1-7 3-1.png');
    % save('H:\MRI\论文材料2\第四章实验\自动分割1135-T1-7 3-1.mat','lc');
    % figure,imshow(imzhuiti,[]);%分割椎体
    % imwrite(imzhuiti,'G:\MRI\第三次组会汇报\实验二\our\k = 100 自动分割 金\自动分割761-T1-7 3-1.png');
    %评价指标
%     bw2 = im2bw(img);%自动分割
    % figure,imshow(bw2);
%     im5 = imread(strcat(mask,'真实',name1));
%     im5 = imread([name2,'.jpg']);%手动分割jpg   png
%     im5 = rgb2gray(im5);
%     figure,imshow(im5);
%     [r,c] = size(im5);
%     bw1 = zeros(r,c);
    % [x,y] = find(img == 255);
    % a = length(x);
%     for i = 1:r
%         for j = 1:c
%             if im5(i,j) >= 247
%                 bw1(i,j) = 1;
%             else
%                 bw1(i,j) = 0;
%             end
%         end
%     end
%     bw1 = im2bw(bw1);
%     figure,imshow(bw1);

%     pos2 = [];%放真实中心点坐标
%     [ls,num]=bwlabel(bw1,4);
%     [ls]=quyupaixu(ls);%只能是6个区域
%     for i = 1:num
%         [xs,ys] = find(ls ==i);
%         pos2(i,1) = round(mean(xs));%纵坐标
%         pos2(i,2) = round(mean(ys));%横坐标
%     end


%     area1=sum(sum(bw1==1));%手动
%     area2=sum(sum(bw2==1));%自动
%     %相减操作
%     imSubtract1=imsubtract(bw1,bw2);%真实
%     imSubtract2=imsubtract(bw2,bw1);%分割
    %转换为二值图像
%     bwSubtract1=im2bw(imSubtract1);
%     bwSubtract2=im2bw(imSubtract2);
%     %求取面积
%     areaSub1=sum(sum(bwSubtract1==1));%真实多余
%     areaSub2=sum(sum(bwSubtract2==1));%分割多余
    %%计算两张图像重叠区域
%     areaOverlap1=area1-areaSub1;%重合
%     areaOverlap2=area2-areaSub2;%重合

%     allarea = areaOverlap1 + areaSub1 + areaSub2;
%     precision = roundn((areaOverlap1 / area1)*100,-2);
%     recall = roundn((areaOverlap1 / area2)*100,-2);
%     Jaccard = roundn((areaOverlap1 / allarea)*100,-2);



    % figure,imshow(imzhuiti,[]);%分割椎体
    hold on;
    %曲度测量
    if op1 == 1
        x11 = pos1(1,2);y11 = pos1(1,1);
        x31 = pos1(3,2);y31 = pos1(3,1);
        x61 = pos1(6,2);y61 = pos1(6,1);
        k1 = (y31-y11)/(x31-x11);
        k2 = (y61-y31)/(x61-x31);
        q1 = qudu(k1,k2);%自动曲度方法1
%         x12 = pos2(1,2);y12 = pos2(1,1);
%         x32 = pos2(3,2);y32 = pos2(3,1);
%         x62 = pos2(6,2);y62 = pos2(6,1);
%         k3 = (y32-y12)/(x32-x12);
%         k4 = (y62-y32)/(x62-x32);
%         q2 = qudu(k3,k4);%手动曲度方法1
%         error1 = abs(q1-q2);
%         bi1 = roundn((error1/q2)*100,-2);
    end
    if op2 == 1
        x13 = pos1(1,2);y13 = pos1(1,1);
        x23 = pos1(2,2);y23 = pos1(2,1);
        x53 = pos1(5,2);y53 = pos1(5,1);
        x63 = pos1(6,2);y63 = pos1(6,1); 
        k5 = (y23-y13)/(x23-x13);
        k6 = (y63-y53)/(x63-x53);
        q3 = qudu(k5,k6);%自动曲度方法2
    %     dian1 =
    %     plot([x13,x13],[y13,y13],'Marker','*','LineWidth',5,'color','b');%画质心
    %     dian2 = plot([x23,x23],[y23,y23],'Marker','*','LineWidth',5,'color','b');
    %     dian3 = plot([x53,x53],[y53,y53],'Marker','*','LineWidth',5,'color','b');
    %     dian4 = plot([x63,x63],[y63,y63],'Marker','*','LineWidth',5,'color','b');


%         x14 = pos2(1,2);y14 = pos2(1,1);
%         x24 = pos2(2,2);y24 = pos2(2,1);
%         x54 = pos2(5,2);y54 = pos2(5,1);
%         x64 = pos2(6,2);y64 = pos2(6,1);
%         k7 = (y24-y14)/(x24-x14);
%         k8 = (y64-y54)/(x64-x54);
%         q4 = qudu(k7,k8);%手动曲度方法2
%         error2 = abs(q3-q4);
%         bi2 = roundn((error2/q4)*100,-2);
    end
    % toc;


%     z = [];
%     z(1,1) = precision;z(1,2) = recall;z(1,3) = Jaccard;z(1,4) = q1;z(1,5) = q2;z(1,6) = error1;
%     z(1,7) = bi1;z(1,8) = q3;z(1,9) = q4;z(1,10) = error2;z(1,11) = bi2;
    mRowRange = num2str(t);
    a = strcat('A',mRowRange);
    z = [];
    z(1,1) = input_age; z(1,2) = input_sex; z(1,3) = q1;z(1,4) = q3;
    xlswrite('曲率.xls',z,'sheet1',a);
    t = t+1;


%     y = [];
%     y(1,1) = m;y(1,2) = seRadius;y(1,3) = pa1;y(1,4) = pa2;y(1,5) = Ec;
% 
% 
%     w = [];
%     w(1,1) = k;w(1,2) = m;w(1,3) = seRadius;w(1,4) = nItr;w(1,5) = Ec;w(1,6) = pa1;w(1,7) = pa2;

end
