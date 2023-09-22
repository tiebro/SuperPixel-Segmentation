function [I] = Pyramid(im)

% %金字塔上采样
% sub11 = im(1:2:end,1:2:end,1);
% sub12 = im(1:2:end,1:2:end,2);
% sub13 = im(1:2:end,1:2:end,3);
% I1=cat(3,sub11,sub12,sub13);
% I1 = lab2rgb(I1);
% figure,imshow(I1);
% 
% sub21 = sub11(1:2:end,1:2:end);
% sub22 = sub12(1:2:end,1:2:end);
% sub23 = sub13(1:2:end,1:2:end);
% I2=cat(3,sub21,sub22,sub23);
% I2 = lab2rgb(I2);
% figure,imshow(I2);
% 
% sub31 = sub21(1:2:end,1:2:end);
% sub32 = sub22(1:2:end,1:2:end);
% sub33 = sub23(1:2:end,1:2:end);
% I3=cat(3,sub31,sub32,sub33);
% I3 = lab2rgb(I3);
% figure,imshow(I3);
% imwrite(I1,'G:\MRI\论文1\金字塔1层23-T1-6 3.jpg');
% imwrite(I2,'G:\MRI\论文1\金字塔2层23-T1-6 3.jpg');
% imwrite(I3,'G:\MRI\论文1\金字塔3层23-T1-6 3.jpg');



imL = im(:,:,1);
ima = im(:,:,2);
imb = im(:,:,3);

%双边滤波
% pa1 = 9;
% pa2 = 30;
% 
% [width, height]=size(imL);
% sigmaSpatial  = min( width, height ) / pa1;%默认16
% samplingSpatial=sigmaSpatial;
% sigmaRange = ( max( imL( : ) ) - min( imL( : ) ) ) / pa2;%默认10
% samplingRange= sigmaRange;
% output = bilateralFilter( imL, imL, sigmaSpatial, sigmaRange, ...
%     samplingSpatial, samplingRange );
% imL = output;
% 
% 
% [width, height]=size(ima);
% sigmaSpatial  = min( width, height ) / pa1;%默认16
% samplingSpatial=sigmaSpatial;
% sigmaRange = ( max( ima( : ) ) - min( imL( : ) ) ) / pa2;%默认10
% samplingRange= sigmaRange;
% output = bilateralFilter( ima, ima, sigmaSpatial, sigmaRange, ...
%     samplingSpatial, samplingRange );
% ima = output;
% 
% 
% [width, height]=size(imb);
% sigmaSpatial  = min( width, height ) / pa1;%默认16
% samplingSpatial=sigmaSpatial;
% sigmaRange = ( max( imb( : ) ) - min( imL( : ) ) ) / pa2;%默认10
% samplingRange= sigmaRange;
% output = bilateralFilter( imb, imb, sigmaSpatial, sigmaRange, ...
%     samplingSpatial, samplingRange );
% imb = output;
% 
% im=cat(3,imL,ima,imb);
% im2 = lab2rgb(im);
% figure,imshow(im2,[]);

% % generate a Gaussian kernel
% ker = fspecial( 'gaussian', 3, 0.5 );
% imL = imfilter( imL, ker, 'conv', 'symmetric', 'same' );
% ima = imfilter( ima, ker, 'conv', 'symmetric', 'same' );
% imb = imfilter( imb, ker, 'conv', 'symmetric', 'same' );
[im_h, im_w] = size( imL );
% image down sampling
imL = imL( 1 : 2 : im_h, 1 : 2 : im_w );
ima = ima( 1 : 2 : im_h, 1 : 2 : im_w );
imb = imb( 1 : 2 : im_h, 1 : 2 : im_w );

I=cat(3,imL,ima,imb);
% im2 = lab2rgb(I);
% figure,imshow(im2,[]);
% imwrite(im2,'G:\MRI\论文1\金字塔1层23-T1-6 3.jpg');