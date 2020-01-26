function [ ssd_gt, rate_gt, RDcost_gt, recDep_gt_hr,ind_lookup ] = mode_gt2( blk, TableGT, rateGT, predictor, QP, BW, bEdge, bSrDep, ind_row,ind_col,recDep )

QStep = QStep_Compute(QP);
lambda = 0.85*2^((QP-6)./3);
ssd = zeros(size(TableGT,3),1);
rate = zeros(size(TableGT,3),1);
RDcost = inf*ones(size(TableGT,3),1);
recDep_lr = zeros(4,4,size(TableGT,3));

for k = 1:size(TableGT,3)
    Evec = TableGT(:,:,k);
    alpha = Evec'*blk(:);    
    alpha_q = floor(abs(alpha)/QStep + 1.0/3).*sign(alpha); 
    vRecDep = Evec*(QStep*alpha_q);
    recDep_lr(:,:,k) = round(reshape(vRecDep,4,4)+predictor(1:2:8,1:2:8));
    ssd(k) = sum(sum((recDep_lr(:,:,k) - bSrDep(1:2:end,1:2:end) ).^2));  % use the ssd of LR block instead for faster implementation

    rateTr = rateGT(k);   % rate of the transform description (the bits to encode the transform index)
    rateTc = AriCod(alpha_q);  % rate of the transform coefficients
    rate(k) = rateTr + rateTc;  % total rate
    RDcost(k) = ssd(k) + lambda.*rate(k);  
end

[minRD,ind] = min(RDcost);
rate_gt = rate(ind);
recDep_gt = recDep_lr(:,:,ind);
recDep_gt_hr = edgeAdaSR_image_2( bEdge,BW,recDep_gt,recDep,ind_row,ind_col,2 );  % bug in this function: the pixels on the last row / column are not well interpolated
ssd_gt = sum(sum((recDep_gt_hr - bSrDep ).^2));
RDcost_gt = ssd_gt + lambda*rate_gt;
ind_lookup = ind;

end

