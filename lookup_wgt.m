function [V, flag, ind] = lookup_wgt(BW,BW_weak,TableGBT)

BW(2:2:end,2:2:end) = 0;
[indX,indY] = find(BW==1);
TableDC = TableGBT(:,1,:);
TableDC(find(abs(TableDC)<0.0001)) = 0;
ind = (1:size(TableGBT,3));

% first check strong edges
if ~isempty(indX)
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

        if isempty(ind_new)
            V = TableGBT(:,:,ind(1));
            ind = ind(1);
            flag = 0;  % find an approximate
%             return;
        else
            ind = ind(ind_new);
        end    
    end
else
    ind = find(sum(abs(TableDC(:,1,:))==0.25)==16);
end

% second, check weak edges
if isempty(indX)
    TableEigen2 = TableGBT(:,2,:);
else
    TableEigen2 = TableGBT(:,3,:);
end
TableEigen2(find(abs(TableEigen2)<0.0001)) = 0;
BW_weak(2:2:end,2:2:end) = 0;
[indX,indY] = find(BW_weak==1);

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
   
    ind_new = find(TableEigen2(uint8(node1),1,ind) .* TableEigen2(uint8(node2),1,ind) < 0);
    
    if isempty(ind_new)
        flag = 0;  % find an approximate
        ind = (1:size(TableGBT,3));
        err = zeros(length(ind),1);
        for iter = 1:length(ind)
            err(iter) = LookupErr(BW,TableDC(:,1,ind(iter)));
        end
        [minErr,ind_approx] = min(err);
        V = TableGBT(:,:,ind(ind_approx));
        ind = ind(ind_approx);            
        return;
        
    else
        ind = ind(ind_new);
    end    
    
%     if isempty(ind_new)
%         V = TableGBT(:,:,ind(1));
%         ind = ind(1);
%         flag = 0;  % find an approximate
%         return;
%     else
%         ind = ind(ind_new);
%     end    
end


% if length(ind) == 1
    V = TableGBT(:,:,ind(1));   
    ind = ind(1);
% else
%     % eliminate unweighted bases
%     for k=1:length(ind)
%         dcNum(k) = length(unique(TableDC(:,1,ind(k))));
%     end
%     [maxNum,ind_uw] = max(dcNum);
%     V = TableGBT(:,:,ind(ind_uw));   
% end
flag = 1;  % find the exactly matched bases

end

