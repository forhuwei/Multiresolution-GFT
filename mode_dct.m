function [ssd,rate,RDcost,alpha_q, recDep] = mode_dct(residue,predictor,bSrDep,QStep,lambda)

% if M == 2
%     residue = residue(1:2:end,1:2:end);
%     predictor = predictor(1:2:end,1:2:end);
% end
alpha = dct2(residue);
alpha_q = floor(abs(alpha)/QStep + 1.0/3).*sign(alpha);
recDep = round(idct2(QStep*alpha_q) + predictor);
% if M == 2
%     recDep = imresize(recDep,2);
% end
        
% RD cost
ssd = sum(sum((recDep-bSrDep).^2));
alpha_zigzag = zigzag(alpha_q);
rateTc = AriCod(alpha_zigzag);
rate = rateTc+1;
% rate_dct = length(find(abs(alpha_zigzag)>=0.0001));
% rate_dct = bSize.^2*entropy_mine(alpha_zigzag);
RDcost = ssd + lambda.*rate;

end