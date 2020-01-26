function [ssd, rate, RDcost, alpha_q, recDep, FlagW,Evec,EdgeW,ind_lookup] = ...
            mode_wgt(residue,bEdge,recDep,TableGBT,predictor,bSrDep,QStep,lambda,BW,ind_row,ind_col,T)
 
    [bSize,bSize] = size(residue);
    lrDep = residue(1:2:bSize,1:2:bSize);

    N1 = 4*(4-1);
    N2 = 4*(4-1);
    N = N1 + N2;
    diff = computeDiff(lrDep);
    alpha = max(diff); 
    rho = 0.1;
    [value,cut] = MRF_hpf(N1,N,4,lrDep,rho,alpha);  % Hochbaum's algorithm 

    weight = 0.13;
    [W,FlagW, EdgeW] = AdjMatW(cut,4,weight,lrDep(:),9);

    if FlagW 
        
        if TableGBT == 0
            ind_lookup = -1;
            D = diag(sum(W,2));
            L = D - W;
            [Evec,Eval] = eig(L);
        else
            [Evec, flag_match, ind_lookup] = lookup_wgt(bEdge,EdgeW,TableGBT);
        end
        
        alpha = Evec'*lrDep(:);       
        alpha_q = floor(abs(alpha)/QStep + 1.0/3).*sign(alpha); 
        vRecDep = Evec*(QStep*alpha_q);
        recDep_lr = round(reshape(vRecDep,4,4)+predictor(1:2:8,1:2:8));
        recDep = edgeAdaSR_image_2( bEdge,BW,recDep_lr,recDep,ind_row,ind_col,2 );

        % RD cost 
        ssd = sum(sum((recDep - bSrDep ).^2));
        rateTr = proxyTrRate(EdgeW);
        rateTc = AriCod(alpha_q);
        rate = rateTr + rateTc;
        RDcost = ssd + lambda.*rate;
    else
        ind_lookup = -1;
        ssd = inf;
        rate = inf;
        RDcost = inf;
        alpha_q = inf;
        recDep = inf;
        Evec = inf;
    end

end