function [ A ] = Cut2AdjMat2(NcutDiscrete1,NcutDiscrete2,indPar1,indPar2)

A = zeros(length(NcutDiscrete1)+length(NcutDiscrete2));

indPart11 = find(NcutDiscrete1 == 1);
for k = 1:length(indPart11)
    A(indPar1(indPart11(k)),indPar1(indPart11(k+1:end))) = 1;
    A(indPar1(indPart11(k+1:end)),indPar1(indPart11(k))) = 1;
end

indPart12 = find(NcutDiscrete1 == 0);
for k = 1:length(indPart12)
    A(indPar1(indPart12(k)),indPar1(indPart12(k+1:end))) = 1;
    A(indPar1(indPart12(k+1:end)),indPar1(indPart12(k))) = 1;
end

indPart21 = find(NcutDiscrete2 == 1);
for k = 1:length(indPart21)
    A(indPar2(indPart21(k)),indPar2(indPart21(k+1:end))) = 1;
    A(indPar2(indPart21(k+1:end)),indPar2(indPart21(k))) = 1;
end

indPart22 = find(NcutDiscrete2 == 0);
for k = 1:length(indPart22)
    A(indPar2(indPart22(k)),indPar2(indPart22(k+1:end))) = 1;
    A(indPar2(indPart22(k+1:end)),indPar2(indPart22(k))) = 1;
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
A = A & MaskA;
A = double(A);

end