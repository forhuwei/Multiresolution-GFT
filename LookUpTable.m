clear all
% clc

%% establish the look-up table by 
% 1. all cases
% load('TableGT/TableGT.mat'); 
% load('TableGT/freq.mat');
load('images\training\TableGT_QP32.mat'); 
load('images\training\freq_QP32.mat'); 

TableGT_all = TableGT;

% 2. sorting
% ind = find(freq>9);  % 94 kinds covering 84%
% ind = find(freq>14);   % QP16: 148 covering 88%
% ind = find(freq>16);   % QP24: 145 covering 88%
ind = find(freq>16);   % QP32: 127 covering 88%

TableGT = TableGT_all(:,:,ind);
freq_s = freq(ind);
length(ind)
sum(freq_s)/sum(freq)
[freq_sorted,IX] = sort(freq_s,'descend');
TableGT_sorted = TableGT(:,:,IX);

len = size(TableGT_sorted,3);
symbols = 1:len;
prob = freq_sorted./sum(freq_sorted);
[dictGT,avglen] = huffmandict(symbols,prob);  % avglen: 2
rateGT = zeros(len,1);
for k = 1:len
    codewords = dictGT{k,2};
    rateGT(k) = length(codewords);
end



%% rearrange the look-up table
% newTable = zeros(size(TableGBT));
% len = size(TableGBT,3);
% M = [1 2; 2 3; 3 4; 5 6; 6 7; 7 8; 9 10; 10 11; 11 12; 13 14; 14 15; 15 16;
%      1 5; 2 6; 3 7; 4 8; 5 9; 6 10; 7 11; 8 12; 9 13; 10 14; 11 15; 12 16];
% 
% ind = 1;
% for iter = 1:24  % 24 edges  
%     for k = 1:len
%         currBases = TableGBT(:,:,k);
%         currDC = currBases(:,1);  % consider only one DC basis
%         if currDC(M(iter,1)) == currDC(M(iter,2))
%             newTable(:,:,ind) = currBases;
%             ind = ind+1;
%         end
%     end   
% end
% flag=1;

%% find the matched or best-approximate GBT
