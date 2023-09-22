% MCLEANUPREGIONS  Morphological clean up of small segments in an image of segmented regions
%                  �ڷָ������ͼ�����γ�С�ε���̬����
% Usage: [seg, Am] = mcleanupregions(seg, seRadius)
%
% Arguments: seg - A region segmented image, such as might be produced by a  
%                  graph cut algorithm.  All pixels in each region are labeled
%                  by an integer.  ����ָ�ͼ�����������ͼ�и��㷨������ÿ�������е��������ض���������ǡ�
%       seRadius - Structuring element radius.  This can be set to 0 in which �ṹԪ�ذ뾶�����������Ϊ0��
%                  case  the function will simply ensure all labeled regions  ����������£��������򵥵�ȷ��
%                  are distinct and relabel them if necessary. ���б�������ǲ�ͬ�ģ����ڱ�Ҫʱ���±�����ǡ�
%
% Returns:   seg - The updated segment image.���µķָ�ͼ��
%             Am - Adjacency matrix of segments.  Am(i, j) indicates whether �ָ���ڽӾ��� Am��i��j��
%                  segments labeled i and j are connected/adjacent     ��ʾ���Ϊi��j�Ķ��Ƿ�����/����
%
% Typical application:����Ӧ��
% If a graph cut or superpixel algorithm fails to converge stray segments���ͼ���и�������㷨�޷������������
% can be left in the result.  This function tries to clean things up by:�ڽ���б�����ɢ�Ρ��˺�������ͨ�����·�ʽ����
% 1) Checking there is only one region for each segment label. If there is���ÿ���α�ǩֻ��һ������
%    more than one region they are given unique labels. ����ж���������Ϊ�����ṩΨһ��ǩ��
% 2) Eliminating regions below the structuring element size �����ṹԪ�ش�С���µ�����
%
% Note that regions labeled 0 are treated as a 'privileged' background region��ע�⣬���Ϊ0
% and is not processed/affected by the function.��������Ϊ����Ȩ���������򣬲��Ҳ��ܸú����Ĵ���/Ӱ�졣
%
% See also: REGIONADJACENCY, RENUMBERREGIONS, CLEANUPREGIONS, MAKEREGIONSDISTINCT

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
%
% March   2013 
% June    2013  Improved morphological cleanup process using distance map

function [seg, Am, mask] = mcleanupregions(seg, seRadius)
option = 2;
    % 1) Ensure every segment is distinct 
    [seg, maxlabel] = makeregionsdistinct(seg);%���ֲ�ͬ�����ģ�maxlabel����ǩ��
    
    % 2) Perform a morphological opening on each segment, subtract the opening��ÿ���߶�
    % from the orignal segment to obtain regions to be reassigned to ��ִ����̬���ڣ�
    % neighbouring segments.  ��ԭʼ�߶��м�ȥ���ڣ��Ի��Ҫ����������߶ε�����
    if seRadius
        se = circularstruct(seRadius);   % Accurate and not noticeably slower
                                         % if radius is small ����뾶��С����׼ȷ�Ҳ������Ա���
                                         %�õ�һ��3*3���߼�����
%       se = strel('disk', seRadius, 4);  % Use approximated disk for speedʹ�ý��ƴ���������ٶ�
        mask = zeros(size(seg));

        if option == 1        
            for l = 1:maxlabel
                b = seg == l;
                mask = mask | (b - imopen(b,se));%������һ��ʹ�����������ù⻬���Ͽ���խ�ļ�Ϻ�����ϸ��ͻ���� ��
            end
            
        else   % Rather than perform a morphological opening on every ���ǰ�˳����ÿ��������������ִ��
               % individual region in sequence the following finds separate ��̬���ڣ������������ҵ�����
               % lists of unconnected regions and performs openings on these. ��δ���������б�����Щ����
               % Typically an image can be covered with only 5 or 6 lists of ��ִ�п��ڡ�ͨ����ͼ����Խ���
               % unconnected regions.  Seems to be about 2X speed of option 5��6��δ���������б��ǡ� 
               % 1. (I was hoping for more...) �ƺ���ѡ��1��2���ٶȡ�����ϣ������......��
            list = finddisconnected(seg);
            
            for n = 1:length(list)
                b = zeros(size(seg));
                for m = 1:length(list{n})
                    b = b | seg == list{n}(m);
                end

                mask = mask | (b - imopen(b,se));
            end
        end
        
        % Compute distance map on inverse of mask������ģ��ת�ľ���ͼ
        [~, idx] = bwdist(~mask);%bwdist�������ڼ���Ԫ��֮��ľ��롣
                                %idx��ʾ�ڸ�Ԫ��������������ķ���Ԫ��λ�ã�������ֵ��
        
        % Assign a label to every pixel in the masked area using the label ofʹ��bwdist����Ĳ��������е���
        % the closest pixel not in the mask as computed by bwdist�����صı�ǩ��Ϊ��ģ�����е�ÿ�����ط����ǩ
        seg(mask) = seg(idx(mask));
    end
    
    % 3) As some regions will have been relabled, possibly broken into several����һЩ���򽫱����·ֽ⣬���ܱ��ֳ�
    % parts, or absorbed into others and no longer exist we ensure all regions�������֣������յ����������Ҳ��ٴ��ڣ�
    % are distinct again, and renumber the regions so that they sequentially����ȷ�����������ٴβ�ͬ�������±������
    % increase from 1.  We also need to reconstruct the adjacency matrix toʹ���Ǵ�1��ʼ�������ӡ����ǻ���Ҫ�ؽ�
    % reflect the changed number of regions and their relabeling.�ڽӾ��� ��ӳ���������ı仯�������±�ǡ�

    seg = makeregionsdistinct(seg);
    [seg, minLabel, maxLabel] = renumberregions(seg);
    Am = regionadjacency(seg);    
    
