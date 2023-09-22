function neighbours = regionQuery(Sp, Am, n, Ec)
    
    lw = 1;
    neighbours = [];
    
    % Get indices of all superpixels connected to superpixel n  获取连接到超像素n的所有超像素的索引
    ind = find(Am(n,:));
    
    for i = ind
        
        dE = cmcdifference([Sp(i).L; Sp(i).a; Sp(i).b],...
                           [Sp(n).L; Sp(n).a; Sp(n).b], lw);
        
        if dE < Ec
            neighbours = [neighbours i];     
        end
    end

