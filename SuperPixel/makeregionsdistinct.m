% MAKEREGIONSDISTINCT Ensures labeled segments are distinctȷ����ǵĶ��ǲ�ͬ��
%
% Usage: [seg, maxlabel] = makeregionsdistinct(seg, connectivity)
%
% Arguments: seg - A region segmented image, such as might be produced by a  ����ָ�ͼ�� 
%                  superpixel or graph cut algorithm.  All pixels in each    ��������ɳ����ػ�ͼ���и�
%                  region are labeled by an integer.           �㷨������ÿ�������е��������ض���������ǡ�
%   connectivity - Optional parameter indicating whether 4 or 8 connectedness
%                  should be used.  Defaults to 4. ��ѡ������ָʾ�Ƿ�Ӧʹ��4��8��ͨ�ԡ� Ĭ��Ϊ4��
%
% Returns:   seg - A labeled image where all segments are distinct.���ͼ���������жζ��ǲ�ͬ�ġ�
%       maxlabel - Maximum segment label number.���α�ǩ�š�
%
% Typical application: A graphcut or superpixel algorithm may terminate in a few  ͼ���и�������㷨����
% cases with multiple regions that have the same label.  This function  �ھ�����ͬ��ǩ�Ķ���������������� 
% identifies these regions and assigns a unique label to them.  ��ֹ���˹��ܿ�ʶ����Щ����Ϊ�����Ψһ��ǩ��
%
% See also: SLIC, CLEANUPREGIONS, RENUMBERREGIONS

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

% June 2013


function [seg, maxlabel] = makeregionsdistinct(seg, connectivity)
    
    if ~exist('connectivity', 'var'), connectivity = 4; end     %������=4
    
    % Ensure every segment is distinct but do not touch segments 
    % with a label of 0 ȷ��ÿ���ζ��ǲ�ͬ�ģ�����Ҫ������ǩΪ0�Ķ�
    labels = unique(seg(:))';%uniqueȡ���ϵĲ��ظ�Ԫ�ع��ɵ�����
    maxlabel = max(labels);
    labels = setdiff(labels,0);  % Remove 0 from the label list�ӱ�ǩ�б���ɾ��0
    
    for l = labels
        [bl,num] = bwlabel(seg==l, connectivity);  %L����seg==l��ÿ����ͨ���������ǩ��num������ͨ����ĸ���
        
        if num > 1  % We have more than one region with the same label�����ж�����������ͬ�ı�ǩ
            for n = 2:num
                maxlabel = maxlabel+1;  % Generate a new label�����±�ǩ
                seg(bl==n) = maxlabel;  % and assign to this segment���������ϸ��
            end
        end
    end

