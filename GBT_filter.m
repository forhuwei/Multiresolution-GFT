function [ lpDep ] = GBT_filter( X, bEdge, M )

X = double(X);
[bSize, bSize] = size(X);
A = AdjMat(bEdge);
D = diag(sum(A,2));  
L = D - A;
[V, Lam] = eig(L);
alpha_norm = V'*X(:);
len = length(alpha_norm);
alpha_norm_hat = zeros(size(alpha_norm));
alpha_norm_hat(1:(len/(M^2))) = alpha_norm(1:(len/(M^2)));

lpDep = round(V*(alpha_norm_hat));
lpDep = vec2mat(lpDep,bSize,bSize);




% [h w] = size(X);
% bSize = 4;     % block size for LR depth
% bCol = floor(w/bSize);     % num of blocks in columns and rows
% bRow = floor(h/bSize);
% 
% recDep = zeros(h, w);
% 
% for Bi = 1:bRow
%     for Bj = 1:bCol
%         
%         %% read depth and edges in a block
%         ind_row = (bSize*(Bi-1)+1) : bSize*Bi;   
%         ind_col = (bSize*(Bj-1)+1) : bSize*Bj;  
%         bSrDep = X(ind_row,ind_col);
% 
%         ind_row_2 = 2*min(ind_row)-1:2*max(ind_row); 
%         ind_col_2 = 2*min(ind_col)-1:2*max(ind_col); 
%         bEdge = BW_big(ind_row_2,ind_col_2);
% %         bEdge(2:2:end,2:2:end) = 0;
%         bEdge(:,end) = 0; bEdge(end,:) = 0;
%         
%         A = AdjMat(bEdge);
%         D = diag(sum(A,2));  
%         L = D - A;
%         [V, Lam] = eig(L);
%         alpha_norm = V'*bSrDep(:);
%         alpha_norm_hat = zeros(size(alpha_norm));
%         alpha_norm_hat(1:16) = alpha_norm(1:16);
% 
%         vRecDep = round(V*(alpha_norm_hat));
%         recDep(ind_row,ind_col) = vec2mat(vRecDep,bSize,bSize);
%         
%     end
% end

end

