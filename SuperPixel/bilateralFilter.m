function output = bilateralFilter( data, edge, sigmaSpatial, sigmaRange, ...
    samplingSpatial, samplingRange )%˫���˲���

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
% Bilateral and Cross-Bilateral Filter using the Bilateral Grid.ʹ��˫�������˫�ߺͿ��˫���˲�����
%
% Bilaterally filters the image 'data' using the edges in the image 'edge'.     ʹ��ͼ��'edge'�еı�Ե
% If 'data' == 'edge', then it the standard bilateral filter.    ˫�߹���ͼ��'����'�����'data'=='edge'��
% Otherwise, it is the 'cross' or 'joint' bilateral filter.         �����Ǳ�׼˫���˲���������
% For convenience, you can also pass in [] for 'edge' for the normal    ����'cross'��'joint'˫���˲����� 
% bilateral filter.        ����Ļ�����Ҳ����Ϊ��ͨ��˫���˲�������[]Ϊ'edge'��
%
% Note that for the cross bilateral filter, data does not need to be        ��ע�⣬���ڽ���˫���˲�����
% defined everywhere.  Undefined values can be set to 'NaN'.  However, edge  ����Ҫ���κεط��������ݡ� 
% *does* need to be defined everywhere.   δ�����ֵ��������Ϊ��NaN���� ���ǣ�edge *ȷʵ*��Ҫ���κεط����塣
%
% data and edge should be of the greyscale, double-precision floating point���ݺͱ�ԵӦ������ͬ��С�ĻҶȣ�
% matrices of the same size (i.e. they should be [ height x width ])˫���ȸ�����󣨼�����Ӧ����[��x��]��
%
% data is the only required argument   data��Ψһ����Ĳ���
%
% edgeMin and edgeMax specifies the min and max values of 'edge' (or 'data'    edgeMin��edgeMaxָ��'edge'��
% for the normal bilateral filter) and is useful when the input is in a   ��Сֵ�����ֵ��������˫���˲���
% range that's not between 0 and 1.  For instance, if you are filtering the  ��'data'����������ķ�Χ����0��1
% L channel of an image that ranges between 0 and 100, set edgeMin to 0 and  ֮��ʱ�����ǳ����á����磬���
% edgeMax to 100.    Ҫ����L ͼ���ͨ����Χ��0��100֮�䣬��edgeMin����Ϊ0����edgeMax����Ϊ100��
% 
% edgeMin defaults to min( edge( : ) ) and edgeMax defaults to max( edge( : ) ).   edgeMinĬ��Ϊmin��edge��:)����
% This is probably *not* what you want, since the input may not span the   edgeMaxĬ��Ϊmax��edge��:)����
% entire range.   �������* not *����Ҫ�ģ���Ϊ������ܲ����Խ������Χ��
%
% sigmaSpatial and sigmaRange specifies the standard deviation of the space
% and range gaussians, respectively.   sigmaSpatial��sigmaRange�ֱ�ָ���ռ�ͷ�Χ��˹�ı�׼ƫ�
% sigmaSpatial defaults to min( width, height ) / 16
% sigmaRange defaults to ( edgeMax - edgeMin ) / 10.
%
% samplingSpatial and samplingRange specifies the amount of downsampling
% used for the approximation.  Higher values use less memory but are also
% less accurate.  The default and recommended values are:
% samplingSpatial��samplingRangeָ�����ڽ��Ƶ��²������� �ϸߵ�ֵʹ�ý��ٵ��ڴ棬
%��Ҳ��̫׼ȷ�� Ĭ��ֵ�ͽ���ֵΪ��
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

% allocate 3D grid ����3D����
downsampledWidth = floor( ( inputWidth - 1 ) / samplingSpatial ) + 1 + 2 * paddingXY;
downsampledHeight = floor( ( inputHeight - 1 ) / samplingSpatial ) + 1 + 2 * paddingXY;
downsampledDepth = floor( edgeDelta / samplingRange ) + 1 + 2 * paddingZ;

gridData = zeros( downsampledHeight, downsampledWidth, downsampledDepth );
gridWeights = zeros( downsampledHeight, downsampledWidth, downsampledDepth );

% compute downsampled indices �����²���ָ��
[ jj, ii ] = meshgrid( 0 : inputWidth - 1, 0 : inputHeight - 1 );

di = round( ii / samplingSpatial ) + paddingXY + 1;
dj = round( jj / samplingSpatial ) + paddingXY + 1;
dz = round( ( edge - edgeMin ) / samplingRange ) + paddingZ + 1;

% perform scatter (there's probably a faster way than this)   ִ�з�ɢ�����ܱ�����죩ͨ����
% normally would do downsampledWeights( di, dj, dk ) = 1, but we have to   ִ��
% perform a summation to do box downsampling    downsampledWeights��di��dj��dk��= 1�������Ǳ���ִ������Խ��п��²���
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

% convolve���
blurredGridData = convn( gridData, kernel, 'same' );
blurredGridWeights = convn( gridWeights, kernel, 'same' );

% divide
blurredGridWeights( blurredGridWeights == 0 ) = -2; % avoid divide by 0, won't read there anyway
normalizedBlurredGrid = blurredGridData ./ blurredGridWeights;
normalizedBlurredGrid( blurredGridWeights < -1 ) = 0; % put 0s where it's undefined
blurredGridWeights( blurredGridWeights < -1 ) = 0; % put zeros back

% upsample�ϲ���
[ jj, ii ] = meshgrid( 0 : inputWidth - 1, 0 : inputHeight - 1 ); % meshgrid does x, then y, so output arguments need to be reversed
                                                                    %meshgrid��x��Ȼ����y���������������Ҫ��ת
% no rounding
di = ( ii / samplingSpatial ) + paddingXY + 1;
dj = ( jj / samplingSpatial ) + paddingXY + 1;
dz = ( edge - edgeMin ) / samplingRange + paddingZ + 1;

% interpn takes rows, then cols, etc       interpnռ���У�Ȼ����cols��
% i.e. size(v,1), then size(v,2), ...
output = interpn( normalizedBlurredGrid, di, dj, dz );