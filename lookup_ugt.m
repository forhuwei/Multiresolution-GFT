function [V, flag, ind_lookup] = lookup_ugt(BW,TableGBT)

BW(2:2:end,2:2:end) = 0;
[indX,indY] = find(BW==1);
TableDC = TableGBT(:,1,:);
TableDC(abs(TableDC)<0.0001) = 0;
ind = (1:size(TableGBT,3));
flag = 1;

for k = 1:length(indX)
    
    if mod(indX(k),2) == 0  % disconnect two vertically adjacent nodes
        nodeX = indX(k)/2;
        nodeY = (indY(k)+1)/2;  
        node1 = 4*(nodeY-1) + nodeX;
        node2 = node1+1;
    else
        nodeX = (indX(k)+1)/2;
        nodeY = indY(k)/2;
        node1 = 4*(nodeY-1) + nodeX;
        node2 = node1 + 4;
    end
   
    ind_new = find(TableDC(uint8(node1),1,ind)~=TableDC(uint8(node2),1,ind));
    
%     if isempty(ind_new)
%         V = TableGBT(:,:,ind(1));
%         ind_lookup = ind(1);
%         flag = 0;  % find an approximate
%         return;
%     else
%         ind = ind(ind_new);
%     end    
    
    if isempty(ind_new)
        flag = 0;  % find an approximate
        ind = (1:size(TableGBT,3));
        err = zeros(length(ind),1);
        for iter = 1:length(ind)
            err(iter) = LookupErr(BW,TableDC(:,1,ind(iter)));
        end
        [minErr,ind_approx] = min(err);
        V = TableGBT(:,:,ind(ind_approx));
        ind_lookup = ind(ind_approx);            
        return;
        
    else
        ind = ind(ind_new);
    end    
end

if length(ind) == 1
    V = TableGBT(:,:,ind);   
    ind_lookup = ind;
else
    % eliminate weighted bases
    for k=1:length(ind)
        dcNum(k) = length(unique(TableDC(:,1,ind(k))));
    end
    [minNum,ind_uw] = min(dcNum);
    V = TableGBT(:,:,ind(ind_uw)); 
    ind_lookup = ind(ind_uw);
end
% flag = 1;  % find the exactly matched bases

end

