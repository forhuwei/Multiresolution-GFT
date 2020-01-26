function [ J ] = GBT_filter_L( X, BW_big, L )

X = double(X);
[h w] = size(X);
wsize = 2*L-1;
J = X;

for i = L:h-L+1
    for j = L:w-L+1
        
%         if i==365 && j==220
%             check=1;
%         end
        
        ind_row = i-(L-1) : i+(L-1);   
        ind_col = j-(L-1) : j+(L-1);  
        X_ij = X(ind_row,ind_col); 
        
        ind_row_2 = 2*min(ind_row)-1:2*max(ind_row); 
        ind_col_2 = 2*min(ind_col)-1:2*max(ind_col); 
        BWb = BW_big(ind_row_2, ind_col_2);  % extract a block of edge map to BWb 
        BWb(2*wsize,:) = 0;   % remove edge pixels along the boundary 
        BWb(:,2*wsize) = 0;
           
        if sum(BWb(:)) ~= 0 
            A_ij = AdjMat(BWb);  
            ind = find(A_ij(wsize*(L-1)+L,:)==1);
            if isempty(ind)
               continue;
            end

            used = zeros(256,1);
            iused = 0;
            dep = 1;
            while(dep<=L-1)
               for k = 1:length(ind)
                   iused = iused + 1;
                   used(iused,1) = ind(k);
                   indexT = find(A_ij(ind(k),:)==1);
                   if isempty(indexT)
                       break;
                   else
                       keep = find(indexT ~= (wsize*(L-1)+L));
                       len = length(keep);
                       m = 1;
                       while(m<=len)
                           index = find(indexT(keep(m)) == used);
                           if index > 0
                              keep(m) = []; 
                              len = len - 1;
                           else
                               m = m+1;
                           end
                       end
                       if dep == 1 && k == 1
                           contain = indexT(keep);
                       else
                           for n=1:length(keep) 
                               contain(1,end+1) = indexT(keep(n));
                          end
                       end
                   end
               end
               if isempty(indexT)
                   break;
               end
               ind = unique(contain);
               dep = dep+1;       
            end

            contain = unique(contain);
            if ~isempty(contain)
                weight = 0.3807*ones(length(contain),1)/(length(contain));
                J(i,j) = 0.6193*X(i,j) + sum(X_ij(contain)*weight);
            end
        else
           H = fspecial('gaussian',[wsize wsize]);
           J(i,j) = sum(sum(X_ij.*H));
        end
        
    end
end

end

