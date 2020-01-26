function [ X_sr ] = edgeAdaSR( BW_big,X_hat,K )
% edge adaptive interpolation and extrapolation w/ HR edges

[imHeight, imWidth] = size(X_hat);
wsize = (2*K-1);
X_sr = 300*ones(K*imHeight,K*imWidth);
X_sr(1:K:K*imHeight, 1:K:K*imWidth) = X_hat;

% BWb(2:2:end,2:2:end) = 0;

for i = 1:K*imHeight
    for j = 1:K*imWidth
        
        if mod(i,K) == 1 && mod(j,K) == 1 
           continue;
        else
            
            ind_row = max(1,i-(K-1)) : min(K*imHeight,i+(K-1));   
            ind_col = max(1,j-(K-1)) : min(K*imWidth,j+(K-1));
            nRow = length(ind_row);
            nCol = length(ind_col);
            X_ij = X_sr(ind_row,ind_col); % interpolation window
            
            ii = K; 
            jj = K;
            if i-(K-1) < 1
                ii = i;
            end
            if j-(K-1) < 1
                jj = j;
            end

            ind_row_2 = 2*min(ind_row)-1:2*max(ind_row); 
            ind_col_2 = 2*min(ind_col)-1:2*max(ind_col); 
            BWb = BW_big(ind_row_2, ind_col_2); 

            % find connected adjacent neighbors within the window
            A_ij = AdjMat(BWb);
%             neighbors = zeros(1,K^2-1);
            col = cat(2,1:(nRow*(jj-1)+ii-1),(nRow*(jj-1)+ii+1):nRow*nCol);
            
            ind = find(A_ij(nRow*(jj-1)+ii,col)~=0); 
            iKeep = find(X_ij(col(ind))~=300);
            len = length(iKeep);
            currInd = col(ind(iKeep));
            neighbors = currInd;
            
            if ~isempty(neighbors)

                X_ij(ii,jj) = 0;
                X_ij(ii,jj) = X_ij(ii,jj)+sum(X_ij(neighbors(1:len)))./len;
                
            else
                ind = find(A_ij(nRow*(jj-1)+ii,:)==0);
                ikeep = find(X_ij(ind) ~= 300);
                [r c] = ind2sub([nRow,nCol],ind(ikeep));
                tmp = (r-ii).*(r-ii)+(c-jj).*(c-jj);
                [tmp IX] = sort(tmp,'ascend');
                X_ij(ii,jj) = X_ij(r(IX(1)),c(IX(1)));
            end
            X_sr(i,j) = X_ij(ii,jj);
        end
       
    end
end

end

