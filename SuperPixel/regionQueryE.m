%------------------------------------------------------------------------
% Find indices of all superpixels adjacent to superpixel n with mean 
% Euclidean colour difference less than Ec.�ҵ��볬����n���ڵ����г����ص�ָ��������ƽ��ŷ�����ɫ��С��Ec��
%
% Arguments:
%             Sp - The struct array of superpixel attributes���������Ե�struct����
%             An - Adjacency matrix�ڽӾ���
%              n - Index of point of interest���ڿ��ǳ����ص�����
%             Ec - Colour distance thresholdɫ�ʾ�����ֵ

function neighbours = regionQueryE(Sp, Am, n, Ec)
    
    E2 = Ec^2;   
    neighbours = [];
    
    % Get indices of all superpixels connected to superpixel n  ��ȡ���ӵ�������n�����г����ص�����
    ind = find(Am(n,:));
    
    for i = ind
        % Test if distance^2 < E^2 
        v = Sp(i).value - Sp(n).value;

        if v'*v < E2 
            neighbours = [neighbours i];     
        end
    end
    
    

