function [ X_sr ] = edgeAdaSR_image( BW_big,X_hat,L )
% edge adaptive interpolation and extrapolation w/ HR edges

% BW_big: 8-connected edge map

[imHeight, imWidth] = size(X_hat);
% wsize = (2*L+1);
X_sr = 300*ones(L*imHeight,L*imWidth);
X_hat = double(X_hat);
X_sr(1:L:L*imHeight, 1:L:L*imWidth) = X_hat;

for i = 1:L*imHeight
    for j = 1:L*imWidth
        
        if mod(i,L) == 1 && mod(j,L) == 1 
           continue;
        else
            
            ind_row = max(1,i-5) : min(L*imHeight,i+5);   
            ind_col = max(1,j-5) : min(L*imWidth,j+5);
            nRow = length(ind_row);
            nCol = length(ind_col);
            X_ij = X_sr(ind_row,ind_col); % interpolation window
            
            ii = 6; 
            jj = 6;
            if i-5 < 1
                ii = i;
            end
            if j-5 < 1
                jj = j;
            end
            ind_row_2 = 2*min(ind_row)-1:2*max(ind_row); 
            ind_col_2 = 2*min(ind_col)-1:2*max(ind_col); 
%             ind_row_2 = 2*min(ind_row)-1:2*ii-1; 
%             ind_col_2 = 2*min(ind_col)-1:2*jj-1; 
            BWb = BW_big(ind_row_2, ind_col_2); 
            BWb(:,end) = 0; BWb(end,:) = 0;
            
            if sum(BWb(:)) == 0
                ind = find(X_ij~=300);
                len = length(ind);
                X_ij(ii,jj) = sum(X_ij(ind(:)))./len;
            
            else
                % find connected adjacent neighbors within the window
                A = AdjMat(BWb);
                ind = find(A(nRow*(jj-1)+ii,:)~=0); 
                if ~isempty(ind)
                    depth = 1;
                    A1 = A;
                    while(depth<=6)
                        ikeep = find(X_ij(ind)~=300);
                        if ~isempty(ikeep)
                            neighbors = ind(ikeep);
                            len = length(neighbors);
                            X_ij(ii,jj) = sum(X_ij(neighbors(1:len)))./len;
                            break;
                        else
                            A1 = A1*A;
                            depth = depth+1;
                            ind = find(A1(nRow*(jj-1)+ii,:)~=0); 
                        end
                    end
                    if depth > 6
                        ind = find(A(nRow*(jj-1)+ii,:)==0);
                        ikeep = find(X_ij(ind) ~= 300);
                        [r c] = ind2sub([nRow,nCol],ind(ikeep));
                        tmp = (r-ii).*(r-ii)+(c-jj).*(c-jj);
                        [tmp IX] = sort(tmp,'ascend');
                        X_ij(ii,jj) = X_ij(r(IX(1)),c(IX(1)));
                    end
                else  % isolated
                    ind = find(A(nRow*(jj-1)+ii,:)==0);
                    ikeep = find(X_ij(ind) ~= 300);
                    [r c] = ind2sub([nRow,nCol],ind(ikeep));
                    tmp = (r-ii).*(r-ii)+(c-jj).*(c-jj);
                    [tmp IX] = sort(tmp,'ascend');
                    X_ij(ii,jj) = X_ij(r(IX(1)),c(IX(1)));
                end
            end
            X_sr(i,j) = X_ij(ii,jj);
        end
    end
end

end

