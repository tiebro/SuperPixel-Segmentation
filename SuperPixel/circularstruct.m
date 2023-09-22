% CIRCULARSTRUCT
%
% Function to construct a circular structuring element
% for morphological operations.  ����������̬ѧ�����Բ�νṹԪ�صĹ��ܡ�
%
% function strel = circularstruct(radius)
%
% Note radius can be a floating point value though the resulting
% circle will be a discrete approximation ע��뾶�����Ǹ���ֵ�������Բ������ɢ����ֵ
%
% Peter Kovesi   March 2000

function strel = circularstruct(radius)

if radius < 1
  error('radius must be >= 1');
end

dia = ceil(2*radius);  % Diameter of structuring element�ṹԪ�ص�ֱ��
                        %ceil�����������ȡ��

if mod(dia,2) == 0     % If diameter is a odd value���ֱ��������ֵ
 dia = dia + 1;        % add 1 to generate a `centre pixel'��1�����ɡ��������ء�
end

r = fix(dia/2);%���㷽��ȡ��
[x,y] = meshgrid(-r:r);
rad = sqrt(x.^2 + y.^2);  %x.^2ÿ��Ԫ��ƽ��
strel = rad <= radius;

