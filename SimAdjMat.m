function [ W, scale_sig ] = SimAdjMat( data, scale_sig, order )

% [BW_big, T, gd] = edge_double_grid_alternate_tree_gd(lrDep,0); 
% alpha = 0.2;  
% sigma = alpha*(0 + max(max(abs(gd))));
% SimMat = exp(-gd.^2/(sigma^2));
% W = AdjMatWei(SimMat);

data = data(:);
distances = zeros(length(data),length(data));
% for j = 1:length(data),
%   distances(j,:) = (sqrt((data(1,:)-data(1,j)).^2 +...
%                 (data(2,:)-data(2,j)).^2));
% end
for i = 1:length(data),
    for j = 1:length(data),
        distances(i,j) = (sqrt((data(i)-data(j)).^2 ));
    end
end

% distances = X2distances(data');

if (~exist('scale_sig')),
    scale_sig = 0.05*max(distances(:)); 
end

if (~exist('order')),
  order = 2;
end

tmp = (distances/scale_sig).^order;

W = exp(-tmp);

end

