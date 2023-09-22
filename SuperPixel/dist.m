%-- dist -------------------------------------------
%
% Usage:  D = dist(C, im, r1, c1, S, m)
% 
% Arguments:   C - Cluster being considered���ڿ��Ǽ�Ⱥ
%             im - sub-image surrounding cluster centreΧ�Ƽ�Ⱥ���ĵ���ͼ��
%         r1, c1 - row and column of top left corner of sub image within the
%                  overall image.   ����ͼ������ͼ�����Ͻǵ��к��С�
%              S - grid spacing ������
%              m - weighting factor between colour and spatial differences.��ɫ�Ϳռ����֮��ļ�Ȩ���ӡ�
%
% Returns:     D - Distance image giving distance of every pixel in the
%                  subimage from the cluster centre  ����ͼ�������ͼ����ÿ�����ؾ���������ĵľ���
%
% Distance = sqrt( dc^2 + (ds/S)^2*m^2 )
% where:
% dc = sqrt(dl^2 + da^2 + db^2)  % Colour distance
% ds = sqrt(dx^2 + dy^2)         % Spatial distance
%
% m is a weighting factor representing the nominal maximum colour distancem�Ǳ�ʾԤ�ڵı�������ɫ����ļ�Ȩ
% expected so that one can rank colour similarity relative to distance���ӣ���˿�������ھ��������Զ���ɫ������
% similarity.  try m in the range [1-40] for L*a*b* space�������� ����L* a* b*�ռ䣬������[1-40]��Χ�ڵ�m
%
% ?? Might be worth trying the Geometric Mean instead ?? ����ֵ�ó��Լ���ƽ��ֵ
%  Distance = sqrt(dc * ds)
% but having a factor 'm' to play with is probably handy  ����һ������'m'���ܷܺ���

% This code could be more efficient  �˴�����ܸ���Ч

function D = dist(C, im, r1, c1, S, m)

    % Squared spatial distance  ƽ���Ŀռ����
    %    ds is a fixed 'image' we should be able to exploit this  ds��һ���̶���'ͼ��'������Ӧ���ܹ�
    %    and use a fixed meshgrid for much of the time somehow...  ����������ĳ�̶ֳ���ʹ�ù̶�������...
    [rows, cols, chan] = size(im);
    [x,y] = meshgrid(c1:(c1+cols-1), r1:(r1+rows-1));
    x = x-C(4);  % x and y dist from cluster centre   x��yԶ�뼯Ⱥ����
    y = y-C(5);
    ds2 = x.^2 + y.^2;
    
    % Squared colour differenceƽ��ɫ��
    for n = 1:3
        im(:,:,n) = (im(:,:,n)-C(n)).^2;
    end
    dc2 = sum(im,3);
    
    D = sqrt(dc2 + ds2/S^2*m^2);

