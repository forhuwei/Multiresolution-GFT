function [ X_sr ] = edgeAdaSR_image_2( bEdge,BW_big,X_hat,recDep,rows,cols,L )
% edge adaptive interpolation and extrapolation w/ HR edges

% BW_big: 8-connected edge map
[imHeight, imWidth] = size(recDep);
[h, w] = size(X_hat);
% wsize = (2*L+1);
X_sr = 300*ones(L*h,L*w);
X_hat = double(X_hat);
X_sr(1:L:L*h, 1:L:L*w) = X_hat;
recDep(rows,cols) = X_sr;

if sum(bEdge(:)) == 0
    X_sr = imresize(X_hat,L);
else
    for i = rows(1):rows(end)
        for j = cols(1):cols(end)
            if mod(i-1,L) == 0 && mod(j-1,L) == 0
                continue;
            end
%             ind_row = max(1,i-3) : min(imHeight,i+3);   
%             ind_col = max(1,j-3) : min(imWidth,j+3);
            ind_row = max(1,i-5) : min(imHeight,i+5);   
            ind_col = max(1,j-5) : min(imWidth,j+5);
            nRow = length(ind_row);
            nCol = length(ind_col);
            X_ij = recDep(ind_row,ind_col); % interpolation window
            
%             ii = 4; 
%             jj = 4;
%             if i-3 < 1
%                 ii = i;
%             end
%             if j-3 < 1
%                 jj = j;
%             end
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
            BWb = BW_big(ind_row_2, ind_col_2); 
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
                        recDep(i,j) = 0;
                        recDep(i,j) = recDep(i,j)+sum(X_ij(neighbors(1:len)))./len;
                        X_ij(ii,jj) = recDep(i,j);
                        break;
                    else
                        A1 = A1*A;
                        depth = depth+1;
                        ind = find(A1(nRow*(jj-1)+ii,:)~=0); 
                    end
                end
            else  % isolated
                ind = find(A(nRow*(jj-1)+ii,:)==0);
                ikeep = find(X_ij(ind) ~= 300);
                [r c] = ind2sub([nRow,nCol],ind(ikeep));
                tmp = (r-ii).*(r-ii)+(c-jj).*(c-jj);
                [tmp IX] = sort(tmp,'ascend');
                recDep(i,j) = X_ij(r(IX(1)),c(IX(1)));
                if recDep(i,j) == 300
                    check = 1;
                end
            end
            
        end
    end
    X_sr = recDep(rows,cols);
end
    
X_sr = round(X_sr);    


end

