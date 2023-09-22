% Find indices of all superpixels adjacent to superpixel n with mean colour
% difference less than Ec.  Use CMC colour difference measure
%   �ҵ��볬����n���ڵ����г����ص���������ƽ��ɫ��С��Ec�� ʹ��CMCɫ�����
% Arguments:
%             Sp - The struct array of superpixel attributes ���������Ե�struct����
%             An - Adjacency matrix �ڽӾ���
%              n - Index of superpixel being considered ���ڿ��ǳ����ص�����
%             Ec - Colour distance threshold ɫ�ʾ�����ֵ

function neighbours = regionQueryCMC(Sp, Am, n, Ec)
    
    lw = 1;
    neighbours = [];
    
    % Get indices of all superpixels connected to superpixel n  ��ȡ���ӵ�������n�����г����ص�����
    ind = find(Am(n,:));
    
    for i = ind

        dE = cmcdifference([Sp(i).L; Sp(i).a; Sp(i).b],...
                           [Sp(n).L; Sp(n).a; Sp(n).b], lw);
        
        if dE < Ec
            neighbours = [neighbours i];     
        end
    end