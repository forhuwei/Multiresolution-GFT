function [V, flag, ind_lookup] = lookup(BW,TableGBT)

BW(2:2:end,2:2:end) = 0;
[indX,indY] = find(BW==1);
TableDC = TableGBT(:,1:3,:);
TableDC(abs(TableDC)<0.0001) = 0;
ind = (1:size(TableGBT,3));

err = zeros(length(ind)-1,1);
for iter = 2:length(ind)
    err(iter-1) = LookupErr(BW,TableDC(:,1:3,iter)); % BW --> Wij
    if err(iter-1)==0       
        ind = iter;
        V = TableGBT(:,:,ind);
        ind_lookup = ind;
        flag = 1;
        return;
    end
end

[minErr,ind_approx] = min(err);
V = TableGBT(:,:,ind_approx+1);
ind_lookup = ind_approx+1; 
flag = 0;

% err = LookupErr(BW_basis,s(:));

end

