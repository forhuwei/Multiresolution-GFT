function [ BW ] = cut2edge(Ncut)

bSize = sqrt(length(Ncut));
blk = reshape(Ncut,bSize,bSize);
[BW,T] = edge_double_grid_alternate_tree(blk,1);

end