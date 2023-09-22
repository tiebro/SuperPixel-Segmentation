clear;clc;close all;
name='图片5';
img=imread([name,'.png']);
figure,imshow(img);


% im = double(img);%读图
% im_R=im(:,:,1);
% im_G= im(:,:,2);
% im_B= im(:,:,3);
% f = [0 -1 0; -1 4 -1; 0 -1 0]/5;
% z_R = filter2(f, im_R);%对原图做高通滤波
% z_G = filter2(f, im_G);%对原图做高通滤波
% z_B = filter2(f, im_B);%对原图做高通滤波
% z = cat(3,z_R,z_G,z_B);
% figure,imshow(z);
% 
% lambda = 0.5;
% result = im2uint8(im + lambda .* z);%将滤波结果与原图相加。
% figure,imshow(result);


%边缘增强
% k=15;           %导热系数,控制平滑
% lambda=0.15;    %控制平滑
% N=20;           %迭代次数
% img=double(img);
% imshow(img,[]);
% img_R=img(:,:,1);
% img_G= img(:,:,2);
% img_B= img(:,:,3);
% [m n]=size(img_R);
% 
% patch=1;
% xxr = zeros(m+patch*2,n+patch*2);
% xxg = zeros(m+patch*2,n+patch*2);
% xxb = zeros(m+patch*2,n+patch*2);
% xxr((patch+1):(m+patch),(patch+1):(n+patch)) = img_R(:,:);
% xxg((patch+1):(m+patch),(patch+1):(n+patch)) = img_G(:,:);
% xxb((patch+1):(m+patch),(patch+1):(n+patch)) = img_B(:,:);
% 
% imgn=zeros(m,n);
% for i=1:N
% 
%     for p=2:m
%         for q=2:n
%             %当前像素的散度，对四个方向分别求偏导，局部不同方向上的变化量，
%             %如果变化较多，就证明是边界，想方法保留边界
%             NI=xxr(p,q)-xxr(p-1,q);
%             SI=xxr(p,q)-xxr(p+1,q);
%             EI=xxr(p,q)-xxr(p,q-1);
%             WI=xxr(p,q)-xxr(p,q+1);
%             
%             %四个方向上的导热系数，该方向变化越大，求得的值越小，从而达到保留边界的目的
%             cN=exp(-NI^2/(k*k));
%             cS=exp(-SI^2/(k*k));
%             cE=exp(-EI^2/(k*k));
%             cW=exp(-WI^2/(k*k));
%             
%             imgn_R(p-1,q-1)=img_R(p-1,q-1)+lambda*(cN*NI+cS*SI+cE*EI+cW*WI);  %扩散后的新值      
%         end
%     end
%     
%     img_R=imgn_R;       %整个图像扩散完毕，用已扩散图像的重新扩散。
% end
% 
% for i=1:N
% 
%     for p=2:m
%         for q=2:n
%             %当前像素的散度，对四个方向分别求偏导，局部不同方向上的变化量，
%             %如果变化较多，就证明是边界，想方法保留边界
%             NI=xxg(p,q)-xxg(p-1,q);
%             SI=xxg(p,q)-xxg(p+1,q);
%             EI=xxg(p,q)-xxg(p,q-1);
%             WI=xxg(p,q)-xxg(p,q+1);
%             
%             %四个方向上的导热系数，该方向变化越大，求得的值越小，从而达到保留边界的目的
%             cN=exp(-NI^2/(k*k));
%             cS=exp(-SI^2/(k*k));
%             cE=exp(-EI^2/(k*k));
%             cW=exp(-WI^2/(k*k));
%             
%             imgn_G(p-1,q-1)=img_G(p-1,q-1)+lambda*(cN*NI+cS*SI+cE*EI+cW*WI);  %扩散后的新值      
%         end
%     end
%     
%     img_G=imgn_G;       %整个图像扩散完毕，用已扩散图像的重新扩散。
% end
% 
% for i=1:N
% 
%     for p=2:m
%         for q=2:n
%             %当前像素的散度，对四个方向分别求偏导，局部不同方向上的变化量，
%             %如果变化较多，就证明是边界，想方法保留边界
%             NI=xxb(p,q)-xxb(p-1,q);
%             SI=xxb(p,q)-xxb(p+1,q);
%             EI=xxb(p,q)-xxb(p,q-1);
%             WI=xxb(p,q)-xxb(p,q+1);
%             
%             %四个方向上的导热系数，该方向变化越大，求得的值越小，从而达到保留边界的目的
%             cN=exp(-NI^2/(k*k));
%             cS=exp(-SI^2/(k*k));
%             cE=exp(-EI^2/(k*k));
%             cW=exp(-WI^2/(k*k));
%             
%             imgn_B(p-1,q-1)=img_B(p-1,q-1)+lambda*(cN*NI+cS*SI+cE*EI+cW*WI);  %扩散后的新值      
%         end
%     end
%     
%     img_B=imgn_B;       %整个图像扩散完毕，用已扩散图像的重新扩散。
% end
% 
% imgn = cat(3,img_R,img_G,img_B);
% figure;
% imshow(imgn,[]);



img = rgb2gray(img);

k=15;           %导热系数,控制平滑
lambda=0.15;    %控制平滑
N=20;           %迭代次数
img=double(img);
imshow(img,[]);
[m n]=size(img);
 
imgn=zeros(m,n);
for i=1:N
 
    for p=2:m-1
        for q=2:n-1
            %当前像素的散度，对四个方向分别求偏导，局部不同方向上的变化量，
            %如果变化较多，就证明是边界，想方法保留边界
            NI=img(p-1,q)-img(p,q);
            SI=img(p+1,q)-img(p,q);
            EI=img(p,q-1)-img(p,q);
            WI=img(p,q+1)-img(p,q);
            
            %四个方向上的导热系数，该方向变化越大，求得的值越小，从而达到保留边界的目的
            cN=exp(-NI^2/(k*k));
            cS=exp(-SI^2/(k*k));
            cE=exp(-EI^2/(k*k));
            cW=exp(-WI^2/(k*k));
            
            imgn(p,q)=img(p,q)+lambda*(cN*NI+cS*SI+cE*EI+cW*WI);  %扩散后的新值      
        end
    end
    
    img=imgn;       %整个图像扩散完毕，用已扩散图像的重新扩散。
end
 
figure;
imshow(imgn,[]);