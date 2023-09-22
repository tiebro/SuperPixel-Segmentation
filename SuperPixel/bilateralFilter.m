function output = bilateralFilter( data, edge, sigmaSpatial, sigmaRange, ...
    samplingSpatial, samplingRange )%双边滤波器

% original src: http://people.csail.mit.edu/jiawen/software/bilateralFilter.m
% original author: Jiawen (Kevin) Chen 
%   <jiawen@csail.mit.edu>
%   http://people.csail.mit.edu/jiawen/
%
 
% output = bilateralFilter( data, edge, ...
%                          edgeMin, edgeMax, ...
%                          sigmaSpatial, sigmaRange, ...
%                          samplingSpatial, samplingRange )
%
% Bilateral and Cross-Bilateral Filter using the Bilateral Grid.使用双边网格的双边和跨边双边滤波器。
%
% Bilaterally filters the image 'data' using the edges in the image 'edge'.     使用图像'edge'中的边缘
% If 'data' == 'edge', then it the standard bilateral filter.    双边过滤图像'数据'。如果'data'=='edge'，
% Otherwise, it is the 'cross' or 'joint' bilateral filter.         则它是标准双边滤波器。否则，
% For convenience, you can also pass in [] for 'edge' for the normal    它是'cross'或'joint'双边滤波器。 
% bilateral filter.        方便的话，你也可以为普通的双边滤波器传递[]为'edge'。
%
% Note that for the cross bilateral filter, data does not need to be        请注意，对于交叉双边滤波器，
% defined everywhere.  Undefined values can be set to 'NaN'.  However, edge  不需要在任何地方定义数据。 
% *does* need to be defined everywhere.   未定义的值可以设置为“NaN”。 但是，edge *确实*需要在任何地方定义。
%
% data and edge should be of the greyscale, double-precision floating point数据和边缘应该是相同大小的灰度，
% matrices of the same size (i.e. they should be [ height x width ])双精度浮点矩阵（即它们应该是[高x宽]）
%
% data is the only required argument   data是唯一必需的参数
%
% edgeMin and edgeMax specifies the min and max values of 'edge' (or 'data'    edgeMin和edgeMax指定'edge'的
% for the normal bilateral filter) and is useful when the input is in a   最小值和最大值（或正常双边滤波器
% range that's not between 0 and 1.  For instance, if you are filtering the  的'data'），当输入的范围不在0和1
% L channel of an image that ranges between 0 and 100, set edgeMin to 0 and  之间时，它非常有用。例如，如果
% edgeMax to 100.    要过滤L 图像的通道范围在0到100之间，将edgeMin设置为0，将edgeMax设置为100。
% 
% edgeMin defaults to min( edge( : ) ) and edgeMax defaults to max( edge( : ) ).   edgeMin默认为min（edge（:)），
% This is probably *not* what you want, since the input may not span the   edgeMax默认为max（edge（:)）。
% entire range.   这可能是* not *你想要的，因为输入可能不会跨越整个范围。
%
% sigmaSpatial and sigmaRange specifies the standard deviation of the space
% and range gaussians, respectively.   sigmaSpatial和sigmaRange分别指定空间和范围高斯的标准偏差。
% sigmaSpatial defaults to min( width, height ) / 16
% sigmaRange defaults to ( edgeMax - edgeMin ) / 10.
%
% samplingSpatial and samplingRange specifies the amount of downsampling
% used for the approximation.  Higher values use less memory but are also
% less accurate.  The default and recommended values are:
% samplingSpatial和samplingRange指定用于近似的下采样量。 较高的值使用较少的内存，
%但也不太准确。 默认值和建议值为：
% samplingSpatial = sigmaSpatial
% samplingRange = sigmaRange

if ~exist( 'edge', 'var' )
    edge = data;
end

inputHeight = size( data, 1 );
inputWidth = size( data, 2 );

if ~exist( 'sigmaSpatial', 'var' )
    sigmaSpatial = min( inputWidth, inputHeight ) / 16;
end

edgeMin = min( edge( : ) );
edgeMax = max( edge( : ) );
edgeDelta = edgeMax - edgeMin;

if ~exist( 'sigmaRange', 'var' )
    sigmaRange = 0.1 * edgeDelta;
end

if ~exist( 'samplingSpatial', 'var' )
    samplingSpatial = sigmaSpatial;
end

if ~exist( 'samplingRange', 'var' )
    samplingRange = sigmaRange;
end

if size( data ) ~= size( edge )
    error( 'data and edge must be of the same size' );
end

% parameters
derivedSigmaSpatial = sigmaSpatial / samplingSpatial;
derivedSigmaRange = sigmaRange / samplingRange;

paddingXY = floor( 2 * derivedSigmaSpatial ) + 1;
paddingZ = floor( 2 * derivedSigmaRange ) + 1;

% allocate 3D grid 分配3D网格
downsampledWidth = floor( ( inputWidth - 1 ) / samplingSpatial ) + 1 + 2 * paddingXY;
downsampledHeight = floor( ( inputHeight - 1 ) / samplingSpatial ) + 1 + 2 * paddingXY;
downsampledDepth = floor( edgeDelta / samplingRange ) + 1 + 2 * paddingZ;

gridData = zeros( downsampledHeight, downsampledWidth, downsampledDepth );
gridWeights = zeros( downsampledHeight, downsampledWidth, downsampledDepth );

% compute downsampled indices 计算下采样指数
[ jj, ii ] = meshgrid( 0 : inputWidth - 1, 0 : inputHeight - 1 );

di = round( ii / samplingSpatial ) + paddingXY + 1;
dj = round( jj / samplingSpatial ) + paddingXY + 1;
dz = round( ( edge - edgeMin ) / samplingRange ) + paddingZ + 1;

% perform scatter (there's probably a faster way than this)   执行分散（可能比这更快）通常会
% normally would do downsampledWeights( di, dj, dk ) = 1, but we have to   执行
% perform a summation to do box downsampling    downsampledWeights（di，dj，dk）= 1，但我们必须执行求和以进行框下采样
for k = 1 : numel( dz )
       
    dataZ = data( k ); % traverses the image column wise, same as di( k )
    if ~isnan( dataZ  )
        
        dik = di( k );
        djk = dj( k );
        dzk = dz( k );

        gridData( dik, djk, dzk ) = gridData( dik, djk, dzk ) + dataZ;
        gridWeights( dik, djk, dzk ) = gridWeights( dik, djk, dzk ) + 1;
        
    end
end

% make gaussian kernel
kernelWidth = 2 * derivedSigmaSpatial + 1;
kernelHeight = kernelWidth;
kernelDepth = 2 * derivedSigmaRange + 1;

halfKernelWidth = floor( kernelWidth / 2 );
halfKernelHeight = floor( kernelHeight / 2 );
halfKernelDepth = floor( kernelDepth / 2 );

[gridX, gridY, gridZ] = meshgrid( 0 : kernelWidth - 1, 0 : kernelHeight - 1, 0 : kernelDepth - 1 );
gridX = gridX - halfKernelWidth;
gridY = gridY - halfKernelHeight;
gridZ = gridZ - halfKernelDepth;
gridRSquared = ( gridX .* gridX + gridY .* gridY ) / ( derivedSigmaSpatial * derivedSigmaSpatial ) + ( gridZ .* gridZ ) / ( derivedSigmaRange * derivedSigmaRange );
kernel = exp( -0.5 * gridRSquared );

% convolve卷积
blurredGridData = convn( gridData, kernel, 'same' );
blurredGridWeights = convn( gridWeights, kernel, 'same' );

% divide
blurredGridWeights( blurredGridWeights == 0 ) = -2; % avoid divide by 0, won't read there anyway
normalizedBlurredGrid = blurredGridData ./ blurredGridWeights;
normalizedBlurredGrid( blurredGridWeights < -1 ) = 0; % put 0s where it's undefined
blurredGridWeights( blurredGridWeights < -1 ) = 0; % put zeros back

% upsample上采样
[ jj, ii ] = meshgrid( 0 : inputWidth - 1, 0 : inputHeight - 1 ); % meshgrid does x, then y, so output arguments need to be reversed
                                                                    %meshgrid做x，然后是y，所以输出参数需要反转
% no rounding
di = ( ii / samplingSpatial ) + paddingXY + 1;
dj = ( jj / samplingSpatial ) + paddingXY + 1;
dz = ( edge - edgeMin ) / samplingRange + paddingZ + 1;

% interpn takes rows, then cols, etc       interpn占用行，然后是cols等
% i.e. size(v,1), then size(v,2), ...
output = interpn( normalizedBlurredGrid, di, dj, dz );