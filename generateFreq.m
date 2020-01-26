load('images\training\TableGT_sorted_88.mat');

freq = zeros(1);
newTable = zeros(16,16,1);
for k = 1:size(TableGT,3)
    if mod(k,10) == 0
        k
    end
    bases_ref = single(TableGT(:,:,k));
    if newTable(:,:,1) == zeros(16,16)
        newTable = bases_ref;
        freq = freq + 1;
        continue;
    end
    
    flag = 0;
    for ind = 1:size(newTable,3)
        bases = single(newTable(:,:,ind));
        [newTable,freq,flag_permute] = GBTMatch(newTable,bases_ref,freq);
%         cnt = 0;
%         for col = 1:16
%             if sum(bases(:,col) == bases_ref(:,col))==16 || sum(bases(:,col) == -bases_ref(:,col))==16
%                 cnt = cnt+1;
%             end
%         end
%         if cnt == 16  % merge if the same 
%             freq(ind) = freq(ind) + 1;
%             flag = 1;
%             break;
%         end
        
    end
%     if flag == 0  
%         newTable(:,:,end+1) = bases_ref;
%         freq(end+1) = 1;
%     end
    
end
