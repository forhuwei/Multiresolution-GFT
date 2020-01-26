clear all;
clc;

imdir=dir('images\testing\*.png');
addpath('images\testing\');

% load the GFT table
load('TableGT\TableGT_sorted_freq3.mat');
TableGT = TableGT_sorted;
load('TableGT\rateGT_freq3.mat');

for i=1:length(imdir)

    fprintf('Compressing image %d:',i);
    fprintf(imdir(i).name);
    fprintf('...\n');
    
    origDep = double(imread(imdir(i).name));
    origDep = origDep(:,:,1);
    
    bSize = 8; % block size
    [h, w] = size(origDep);
    bRow = floor(h/bSize);    
    bCol = floor(w/bSize);   
    [BW, T] = edge_double_grid_alternate_tree(origDep,0); % compute strong boundaries in the image
    BW(2:2:end,2:2:end) = 0;
             
    for QP = 16:4:28
        
        QStep = QStep_Compute(QP);   % quantization step size
        lambda = 0.85*2^((QP-6)./3); % the weighting parameter for calculating the RD cost
        MSE = 0; 
        RATE = 0;  
        recDep = origDep;
        
        indTr = zeros(bRow,bCol);        % store transform indices for each block
        
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
                predictor = intra_prediction(origDep,BW,Bi,Bj,bSize);
                residue = bSrDep - predictor;

              %% step 2. construct transforms
                % Transform Mode 1: 8*8 DCT 
                [ssd_dct,rate_dct,RDcost_dct,alpha_dct, recDep_dct] = mode_dct(residue,predictor,bSrDep,QStep,lambda);
                
                % Transform Mode 2: LR-GT (multi-resolution GFT, block downsampled by two)
                LrImg = residue(1:2:end,1:2:end);
                [ ssd_gt, rate_gt, RDcost_gt, recDep_gt,ind_lookup ] = mode_gt2(LrImg,TableGT,rateGT,predictor,QP,BW,bEdge,bSrDep,ind_row,ind_col,recDep);
                 
                % mode decision
                [minRDcost,ind] = min([RDcost_gt RDcost_dct]);

                if ind == 2 
                    MSE = MSE + ssd_dct;
                    RATE = RATE + rate_dct;
                    recDep(ind_row,ind_col) = recDep_dct;
                    
                else
                    MSE = MSE + ssd_gt;
                    RATE = RATE + rate_gt;
                    recDep(ind_row,ind_col) = recDep_gt;
                    indTr(Bi,Bj) = ind_lookup;
                end
               
            end
        end
        
        MSE = MSE/(h*w);
        PSNR = 10*log10(255*255/MSE);
        
        imwrite(uint8(recDep),['results\MR-GFT\',imdir(i).name,'_QP',int2str(QP),'_MRGFT.png']);
        fid=fopen('results\MR-GFT\log.txt','a+');
        fprintf(fid,'%s\n\r',imdir(i).name);
        fprintf(fid,'QP = %f, PSNR = %f, Rate = %f\n\r',QP,PSNR,RATE);
        fclose(fid);      
        save(['results\indTr\',imdir(i).name,'_QP',int2str(QP),'.mat'],'indTr');
        
    end
end






