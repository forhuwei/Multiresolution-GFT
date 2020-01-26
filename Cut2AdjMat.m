function [ A ] = Cut2AdjMat(NcutDiscrete,ind)

A = zeros(length(NcutDiscrete));

indPart1 = find(NcutDiscrete == 1);
for k = 1:length(indPart1)
    A(indPart1(k),indPart1(k+1:end)) = 1;
    A(indPart1(k+1:end),indPart1(k)) = 1;
end

indPart2 = find(NcutDiscrete == 0);
for k = 1:length(indPart2)
    A(indPart2(k),indPart2(k+1:end)) = 1;
    A(indPart2(k+1:end),indPart2(k)) = 1;
end

% covert to 4-connected
MaskA = [0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0;
         1     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0;
         0     1     0     1     0     0     1     0     0     0     0     0     0     0     0     0;
         0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0;
         1     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0;
         0     1     0     0     1     0     1     0     0     1     0     0     0     0     0     0;
         0     0     1     0     0     1     0     1     0     0     1     0     0     0     0     0;
         0     0     0     1     0     0     1     0     0     0     0     1     0     0     0     0;
         0     0     0     0     1     0     0     0     0     1     0     0     1     0     0     0;
         0     0     0     0     0     1     0     0     1     0     1     0     0     1     0     0;
         0     0     0     0     0     0     1     0     0     1     0     1     0     0     1     0;
         0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     1;
         0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0;
         0     0     0     0     0     0     0     0     0     1     0     0     1     0     1     0;
         0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     1;
         0     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0];
MaskA = MaskA(ind,:);
MaskA = MaskA(:,ind);
A = A & MaskA;
A = double(A);


end