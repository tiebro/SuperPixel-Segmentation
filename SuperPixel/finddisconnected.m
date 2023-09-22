% FINDDISCONNECTED find groupings of disconnected labeled regions
%                   �ҵ��Ͽ����ӵı������ķ���
% Usage: list = finddisconnected(l)
%
% Argument:   l - A labeled image segmenting an image into regions, such as  ���ͼ��ͼ��ָ����������
%                 might be produced by a graph cut or superpixel algorithm. ����ͨ��ͼ���и�������㷨������
%                 All pixels in each region are labeled by an integer.      ÿ�������е��������ض���������ǡ�
%
% Returns: list - A cell array of lists of regions that are notδ���ӵ������б�ĵ�Ԫ���С� 
%                 connected. Typically there are 5 to 6 lists.ͨ����5��6���б�
%
% Used by MCLEANUPREGIONS to reduce the number of morphological closing
% operations  ��MCLEANUPREGIONS���ڼ�����̬�պϲ���������
%
% See also: MCLEANUPREGIONS, REGIONADJACENCY

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

% PK July 2013


function list = finddisconnected(l)
 
    debug = 0;
    [Am, Al] = regionadjacency(l);%�����ǵķָ������ͼ����ڽӾ���
    
    N = max(l(:));  % number of labels
    
    % Array for keeping track of visited labels���ڸ��ٷ��ʹ��ı�ǩ������
    visited = zeros(N,1);

    list = {};
    listNo = 0;
    for n = 1:N

        if ~visited(n)
            listNo = listNo + 1;
            list{listNo} = n;
            visited(n) = 1;
            
            % Find all regions not directly connected to n and not visited����δֱ�����ӵ�n��δ���ʵ���������
            notConnected = setdiff(find(~Am(n,:)), find(visited));
            %c = setdiff(A, B)   ������A���У���B��û�е�ֵ��������������������򷵻ء�
            
            % For each unconnected region check that it is not already����ÿ��δ���ӵ����򣬼�����Ƿ���δ
            % connected to a region in the list. If not, add to list���ӵ��б��е����� ���û�У�����ӵ��б���
            for m = notConnected
                if isempty(intersect(Al{m}, list{listNo}))
                    list{listNo} = [list{listNo} m];
                    visited(m) = 1;
                end
            end
         end % if not visited(n)
        
    end
    
    % Display each list of unconncted regions as an image ��ÿ��δ���������б���ʾΪͼ��
    if debug   %����ģʽ
        for n = 1:length(list)
            
            mask = zeros(size(l));
            for m = 1:length(list{n})
                mask = mask | l == list{n}(m);
            end
            
            fprintf('list %d of %d length %d \n', n, length(list), length(list{n}))
            show(mask);
            keypause
        end
    end
    