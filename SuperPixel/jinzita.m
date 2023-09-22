clear;clc;close all;

im = imread( '²Ã¼ô23-T1-6 3.jpg' );
if ndims( im ) == 3 
    im = rgb2gray( im );
end

% generate a Gaussian kernel
ker = fspecial( 'gaussian', 3, 0.5 );
im = imfilter( im, ker, 'conv', 'symmetric', 'same' );
[im_h, im_w] = size( im );

% image down sampling
im2 = im( 1 : 2 : im_h, 1 : 2 : im_w );

figure, imshow( im ); title( 'original image' );
figure, imshow( im2 ); title( 'resampled image' );