function [BW_big, T] = edge_double_grid_alternate_tree(X,flag)

X = double(X);
[Nx,Ny]=size(X);
Y = zeros(2*Nx,2*Ny);

for i = 2:Nx-1
    for j = 2:Ny-1
        Y(2*(i-0.5)-1,2*j-1) = X(i,j)-X(i-1,j);
%         Y(2*(i-0.5)-1,2*(j+0.5)-1) = max(abs(X(i,j)-X(i-1,j+1)),abs(X(i-1,j)-X(i,j+1))); % huwei
        Y(2*i-1,2*(j+0.5)-1) = X(i,j)-X(i,j+1);
%         Y(2*(i+0.5)-1,2*(j+0.5)-1) = max(abs(X(i,j)-X(i+1,j+1)),abs(X(i+1,j)-X(i,j+1))); % huwei
    end
end

% First Column
% for i = 2:Nx-1
for i = 2:Nx

    j = 1;
    
    Y(2*(i-0.5)-1,2*j-1) = X(i,j)-X(i-1,j);
%     Y(2*(i-0.5)-1,2*(j+0.5)-1) = max(abs(X(i,j)-X(i-1,j+1)),abs(X(i-1,j)-X(i,j+1))); % huwei
    Y(2*i-1,2*(j+0.5)-1) = X(i,j)-X(i,j+1);
%     Y(2*(i+0.5)-1,2*(j+0.5)-1) = max(abs(X(i,j)-X(i+1,j+1)),abs(X(i+1,j)-X(i,j+1))); % huwei
end

% Last Column
for i = 2:Nx
    j = Ny;
    
    Y(2*(i-0.5)-1,2*j-1) = X(i,j)-X(i-1,j);
end

% First Row
for j = 1:Ny-1
    i = 1;
    
    Y(2*i-1,2*(j+0.5)-1) = X(i,j)-X(i,j+1);
%     Y(2*(i+0.5)-1,2*(j+0.5)-1) = max(abs(X(i,j)-X(i+1,j+1)),abs(X(i+1,j)-X(i,j+1)));  % huwei
end

% Last Row
for j = 2:Ny-1
    i = Nx;
    
    Y(2*(i-0.5)-1,2*j-1) = X(i,j)-X(i-1,j);
%     Y(2*(i-0.5)-1,2*(j+0.5)-1) = max(abs(X(i,j)-X(i-1,j+1)),abs(X(i-1,j)-X(i,j+1))); % huwei
    Y(2*i-1,2*(j+0.5)-1) = X(i,j)-X(i,j+1);
end

[ifull,jfull] = meshgrid(1:2*Nx,1:2*Ny);
ifull = ifull(:); jfull = jfull(:);
ixkeep = find((mod(ifull,2) == 0) | (mod(jfull,2) == 0));

% iximage = find((mod(ifull,2) == 1) & (mod(jfull,2) == 1));

% X_big = zeros(2*Nx,2*Ny);
% X_big(iximage) = double(X);
% X_big = reshape(X_big,2*Nx,2*Ny);

mu_g = mean(abs(Y(ixkeep)));
std_g = std(abs(Y(ixkeep)));

T = (mu_g + std_g);
% T = 2*(mu_g + std_g);
% T = mu_g+2*std_g;          % threshold for an edge point
if flag == 1
    T = 0;
end
BW_big = zeros(2*Nx, 2*Ny);
BW_big(find(abs(Y)>T)) = 1;
% BW_big(find(abs(Y)>4*T)) = 1;

% [xe,ye] = meshgrid(1:0.5:Nx+0.5,1:0.5:Ny+0.5);
% 
% [xg yg] = meshgrid(2:2:2*Nx, 2:2:2*Ny);

% ix_dual = find((mod(ifull,2) == 0) & (mod(jfull,2) == 0));
% BW_big2=zeros(size(BW_big));
% BW_big2(ix_dual)=BW_big(ix_dual);
% 
% for i=1:2:2*Nx-1
%     v1=find(BW_big(i,:));
%     v2=find(BW_big(:,i));
%     if i>1
%         BW_big2(i-1,v1)=1;
%         BW_big2(v2,i-1)=1;
%     end
%     BW_big2(i+1,v1)=1;
%     BW_big2(v2,i+1)=1;    
% end
% 
% 
% BW_dual = reshape(BW_big2(ix_dual), [Nx Ny]);

% BW2=BW_big(ixkeep);
% BW2 = reshape(BW2,length(BW2)/2,2);
% 
% 
% % [y,nbr_bits] = perform_jbig_coding(BW_big);
% [y,nbr_bits] = perform_jbig_coding(BW_dual);




% [L,num] = bwlabel(BW_big);  
% % [L,num] = bwlabel(BW_big,4);
% 
% %Remove "small" connected components  
% for n = 1:num
%     ixremove = find(L==n);
% 
%     if length(ixremove) <= 10
%         BW_big(ixremove) = 0;
%     end
% end