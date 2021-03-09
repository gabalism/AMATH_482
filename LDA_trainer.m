function [U,S,V,threshold,w,sortones,sorttwos] = LDA_trainer(first_number,second_number,feature)
    
    nd = size(first_number,2);
    nc = size(second_number,2);

    [U,S,V] = svd([first_number second_number],'econ'); 
    numbers = S*V';
    U = U(:,1:feature); % Add this in
    ones = numbers(1:feature,1:nd);
    twos = numbers(1:feature,nd+1:nd+nc);
    md = mean(ones,2);
    mc = mean(twos,2);

    Sw = 0;
    for k=1:nd
        Sw = Sw + (ones(:,k)-md)*(ones(:,k)-md)';
    end
    for k=1:nc
        Sw = Sw + (twos(:,k)-mc)*(twos(:,k)-mc)'; % within class matrix
    end
    Sb = (md-mc)*(md-mc)'; %between class scatter matrix
    
    [V2,D] = eig(Sb,Sw);
    [lambda,ind] = max(abs(diag(D)));
    w = V2(:,ind);
    w = w/norm(w,2);
    vfnumber = w'*ones;
    vsnumber = w'*twos;
    
    if mean(vfnumber)>mean(vsnumber)
        w = -w;
        vfnumber = -vfnumber;
        vsnumber = -vsnumber;
    end
    
    % Don't need plotting here
    sortones = sort(vfnumber);
    sorttwos = sort(vsnumber);
    t1 = length(sortones);
    t2 = 1;
    while sortones(t1)>sorttwos(t2)
    t1 = t1-1;
    t2 = t2+1;
    end
    threshold = (sortones(t1)+sorttwos(t2))/2;

    % We don't need to plot results
end

