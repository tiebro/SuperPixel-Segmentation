function neighbours = regionQuery(Sp, Am, n, Ec)
    
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

