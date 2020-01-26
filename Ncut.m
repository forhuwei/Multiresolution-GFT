function [ NcutDiscrete ] = Ncut(lrDep)

len = length(lrDep);

% 1. compute the weighted adjacency matrix
[W, scale_sig] = SimAdjMat(lrDep);
if scale_sig <= 0.1
    NcutDiscrete = ones(size(lrDep(:)));
else
    % 2. Ncut
    % nbEigenValues = 6;
    % [NcutEigenvectors,NcutEigenvalues] = ncut(W,nbEigenValues);
    D = diag(sum(W,2));
    L = D - W;
    [EigenVectors,EigenValues] = eig(L,D);
    [NcutDiscrete,NcutEigenvectors] = discretisation(EigenVectors(:,1:2));

    % 3. construct GBT
    % to do: decide if further bipartition is necessary?
    NcutDiscrete = full(NcutDiscrete);
    % indBlk = reshape(NcutDiscrete(:,1),bSize,bSize);
    % [BW,T] = edge_double_grid_alternate_tree(indBlk);
    % A = AdjMat(BW);
end

end