function [ssd, rate, RDcost, mRecDep, Evec,ind] = mode_uwgt(residue,predictor,bSrDep,BW,recDep,QP,ind_row,ind_col,TableGBT)

if TableGBT == 0
    ind = -1;
end

[bSize,bSize] = size(residue);
QStep = QStep_Compute(QP);
lambda = 0.85*2^((QP-6)./3);

ind_row_2 = 2*min(ind_row)-1:2*max(ind_row); 
ind_col_2 = 2*min(ind_col)-1:2*max(ind_col); 
bEdge = BW(ind_row_2,ind_col_2);

% BW_s = getBWs( bEdge,M );
% BW_s(2:2:end,2:2:end) = 0; 
% BW_s(:,end) = 0; BW_s(end,:) = 0;

% lpDep = GBT_filter(residue, bEdge, 2);
lrDep = residue(1:2:bSize,1:2:bSize);
% if sum(BW_s(:)) == 0
%     V1 = TableGBT(:,:,1);
% else
%     [V1,flag_match] = lookup(BW_s,TableGBT);
% end

% the first cut
[ NcutDiscrete ] = Ncut(lrDep);
A = Cut2AdjMat(NcutDiscrete(:,1),1:16);
D = diag(sum(A,2));  
L = D - A;
[V1, Lam] = eig(L);

% lookup
% [ edge1 ] = cut2edge(NcutDiscrete(:,1));
% [V1, flag1, ind_lookup1] = lookup_ugt(edge1,TableGBT);


% ind1 = 1;
% ind = 1;
% if length(unique(NcutDiscrete(:,1))) == 1
%     V1 = TableGBT(:,:,1);
%     ind1 = 1;
% else
%     [ BW1 ] = cut2edge(NcutDiscrete(:,1));
%     [V1,flag_match,ind1] = lookup(BW1,TableGBT);
% end
alpha = V1'*lrDep(:); 
alpha_q = floor(abs(alpha)/QStep + 1.0/3).*sign(alpha); 
vRecDep_lr = round(V1*(QStep*alpha_q));
mRecDep_lr = round(reshape(vRecDep_lr,4,4) + predictor(1:2:bSize,1:2:bSize));
mRecDepS = edgeAdaSR_image_2( bEdge,BW,mRecDep_lr,recDep,ind_row,ind_col,2 );
ssd1 = sum(sum((mRecDepS(:) - bSrDep(:) ).^2));
rate_coeff1 = AriCod(alpha_q);
BWb = cut2edge(NcutDiscrete(:,1));
% [rate_coeff1, rateGT_1] = EncodeBlock_CABAC2(alpha_q,BWb,QP);
% rateGT_11 = rateGT(ind_lookup);
% if sum(BWb(:)) ~= 0
%     check = 1;
% end
rateGT_1 = proxyTrRate(BWb);
rate1 = rate_coeff1 + rateGT_1;
RDcostS = ssd1 + lambda.*rate1;

% recursively repartition the segmented parts if necessary
if ~isempty(find(NcutDiscrete(:,1) == 0, 1)) && ~isempty(find(NcutDiscrete(:,1) == 1, 1))
    flag = 1;
    indPar1 = find(NcutDiscrete(:,1) == 1);
    [ NcutDiscrete1 ] = Ncut(lrDep(indPar1));
    indPar2 = find(NcutDiscrete(:,1) == 0);
    [ NcutDiscrete2 ] = Ncut(lrDep(indPar2));
    
    if length(unique(NcutDiscrete1))==1 && length(unique(NcutDiscrete2))==1
        ssd = ssd1;
        rate = rate1;
        RDcost = RDcostS;
        mRecDep = mRecDepS;
        Evec = V1;
%         ind = ind_lookup1;
        return;
    end

    % construct GBT in the entire block
    A = Cut2AdjMat2(NcutDiscrete1(:,1),NcutDiscrete2(:,1),indPar1,indPar2);
    D = diag(sum(A,2));  
    L = D - A;
    [V2, Lam] = eig(L);
    % lookup
    NcutDiscrete2(:,1) = NcutDiscrete2(:,1) + 2;
    NcutDiscrete = zeros(16,1);
    NcutDiscrete(1:length(NcutDiscrete1(:,1))) = NcutDiscrete1(:,1);
    NcutDiscrete(length(NcutDiscrete1(:,1))+1:16) = NcutDiscrete2(:,1); 
    
%     [ edge2 ] = cut2edge(NcutDiscrete);
%     [V2, flag2, ind_lookup2] = lookup_ugt(edge2,TableGBT);
%     ind2 = 1;
%     NcutDiscreteAll = NcutDiscrete(:,1);
%     NcutDiscreteAll(indPar1) = 2*(NcutDiscrete1(:,1)+1);
%     NcutDiscreteAll(indPar2) = 3*NcutDiscrete2(:,1);
%     BW2 = cut2edge(NcutDiscreteAll);
%     [V2,flag_match,ind2] = lookup(BW2,TableGBT);

    alpha = V2'*lrDep(:); 
    alpha_q_M = floor(abs(alpha)/QStep + 1.0/3).*sign(alpha); 
    vRecDep_lr = V2*(QStep*alpha_q_M);
    mRecDep_lr = round(reshape(vRecDep_lr,4,4) + predictor(1:2:bSize,1:2:bSize));
    mRecDepM = edgeAdaSR_image_2( bEdge,BW,mRecDep_lr,recDep,ind_row,ind_col,2 );
    ssd2 = sum(sum((mRecDepM(:) - bSrDep(:) ).^2));
    BWb = cut2edge(NcutDiscrete);
%     [rate_coeff2, rateGT_2] = EncodeBlock_CABAC2(alpha_q,BWb,QP)
    rate_coeff2 = AriCod(alpha_q_M);
    rateGT_2 = proxyTrRate(BWb);
    rate2 = rate_coeff2 + rateGT_2;
    RDcostM = ssd2 + lambda.*rate2;


    % decide if the current partition shuold be subdivided 
    if RDcostM < RDcostS   
        ssd = ssd2;
        rate = rate2;
        RDcost = RDcostM;
        mRecDep = mRecDepM;
        Evec = V2;
        alpha = alpha_q_M;
%         ind = ind2;
    else
        ssd = ssd1;
        rate = rate1;
        RDcost = RDcostS;
        mRecDep = mRecDepS;
        Evec = V1;
        alpha = alpha_q;
%         ind = ind1;
    end
else
    ssd = ssd1;
    rate = rate1;
    RDcost = RDcostS;
    mRecDep = mRecDepS;
    Evec = V1;
    flag = 0;
    alpha = alpha_q;
end

end
