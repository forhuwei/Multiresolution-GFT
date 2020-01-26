clear all;
clc;

imdir=dir('images\training\*.png');
addpath('images\training\');
bSize = 8;

TableGT = zeros(16,16,1);
freq = zeros(1);
    
for i=1:length(imdir)
  
    fprintf('Compressing image %d:',i);
    fprintf(imdir(i).name);
    fprintf('...\n');
    
    origDep = double(imread(imdir(i).name));
    origDep = origDep(:,:,1);
    [h, w] = size(origDep);
    bRow = floor(h/bSize);    
    bCol = floor(w/bSize);   
    [BW, T] = edge_double_grid_alternate_tree(origDep,0);
          
    QP = 16;  
    QStep = QStep_Compute(QP);
    lambda = 0.85*2^((QP-6)./3);

    MSE = 0; 
    RATE = 0;  
    recDep = origDep;

    indTr = zeros(bRow,bCol);        % transform index
    flag_wgt = zeros(bRow,bCol);   % indicate which blocks use wgt
    flag_uwgt = zeros(bRow,bCol);

    for Bi = 1:bRow
        Bi
        for Bj = 1:bCol

            % read image and edges in a block
            ind_row = (bSize*(Bi-1)+1) : bSize*Bi;   
            ind_col = (bSize*(Bj-1)+1) : bSize*Bj;  
            bSrDep = origDep(ind_row,ind_col);

            ind_row_2 = 2*min(ind_row)-1:2*max(ind_row); 
            ind_col_2 = 2*min(ind_col)-1:2*max(ind_col); 
            bEdge = BW(ind_row_2,ind_col_2);
            bEdge(2:2:end,2:2:end) = 0;   
            bEdge(:,end) = 0; bEdge(end,:) = 0;

           %% step 1. intra prediction
            predictor = intra_prediction(recDep,BW,Bi,Bj,bSize);
            residue = bSrDep - predictor;

           %% step 2. construct transforms
            % Mode 1: 8*8 DCT 
            [ssd_dct,rate_dct,RDcost_dct,alpha_dct, recDep_dct] = mode_dct(residue,predictor,bSrDep,QStep,lambda);

            % Mode 2: 4*4 LR-UWGT
            [ssd_uwgt, rate_uwgt, RDcost_uwgt, recDep_uwgt, Evec_uwgt,ind_lookup_uwgt] =...
                mode_uwgt(residue,predictor,bSrDep,BW,recDep,QP,ind_row,ind_col,0);

            % Mode 3: 4*4 LR-WGT
            [ssd_wgt, rate_wgt, RDcost_wgt, alpha_wgt, recDep_wgt,FlagW,Evec_wgt,Edge_wgt,ind_lookup_wgt] = ...
                mode_wgt(residue,bEdge,recDep,0,predictor,bSrDep,QStep,lambda,BW,ind_row,ind_col,T);

            % mode decision
            [minRDcost,ind] = min([RDcost_dct RDcost_uwgt RDcost_wgt]);

            if ind == 1 
                MSE = MSE + ssd_dct;
                RATE = RATE + rate_dct;
                recDep(ind_row,ind_col) = recDep_dct;

            elseif ind == 2
                MSE = MSE + ssd_uwgt;
                RATE = RATE + rate_uwgt;
                recDep(ind_row,ind_col) = recDep_uwgt;
                indTr(Bi,Bj) = ind_lookup_uwgt;
                bases = Evec_uwgt;

            else
                MSE = MSE + ssd_wgt;
                RATE = RATE + rate_wgt;
                recDep(ind_row,ind_col) = recDep_wgt;
                indTr(Bi,Bj) = ind_lookup_wgt;
                bases = Evec_wgt;
            end
            if ind ~= 1
                [TableGT,freq] = GBTMatch(TableGT,bases,freq);
            end
           
        end       
    end
    save('images\training\TableGT_QP16','TableGT');
end






