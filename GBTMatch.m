function [TableGBT,freq,flag_permute] = GBTMatch(TableGBT,bases,freq)
    

flag = 0;
flag_permute = 0;
bases = single(bases);

if TableGBT(:,:,1) == zeros(16,16)
    TableGBT = bases;
    freq = freq + 1;
    return;
end

% ------- permutation matrices -------
P = zeros(16,16,7);
% mirror w.r.t horizontal central line
P(1,4,1) = 1; P(5,8,1) = 1; P(9,12,1) = 1; P(13,16,1) = 1;  
P(4,1,1) = 1; P(8,5,1) = 1; P(12,9,1) = 1; P(16,13,1) = 1;  
P(2,3,1) = 1; P(6,7,1) = 1; P(10,11,1) = 1; P(14,15,1) = 1;
P(3,2,1) = 1; P(7,6,1) = 1; P(11,10,1) = 1; P(15,14,1) = 1;
% mirror w.r.t vertical central line
P(1,13,2) = 1; P(2,14,2) = 1; P(3,15,2) = 1; P(4,16,2) = 1;  
P(13,1,2) = 1; P(14,2,2) = 1; P(15,3,2) = 1; P(16,4,2) = 1;  
P(5,9,2) = 1; P(6,10,2) = 1; P(7,11,2) = 1; P(8,12,2) = 1;
P(9,5,2) = 1; P(10,6,2) = 1; P(11,7,2) = 1; P(12,8,2) = 1;
% mirror w.r.t diagonal central line
P(2,5,3) = 1; P(3,9,3) = 1; P(4,13,3) = 1; P(7,10,3) = 1; P(8,14,3) = 1; P(12,15,3) = 1; 
P(5,2,3) = 1; P(9,3,3) = 1; P(13,4,3) = 1; P(10,7,3) = 1; P(14,8,3) = 1; P(15,12,3) = 1; 
P(1,1,3) = 1; P(6,6,3) = 1; P(11,11,3) = 1; P(16,16,3) = 1;
% mirror w.r.t anti-diagonal central line
P(3,8,4) = 1; P(2,12,4) = 1; P(1,16,4) = 1; P(6,11,4) = 1; P(5,15,4) = 1; P(9,14,4) = 1; 
P(8,3,4) = 1; P(12,2,4) = 1; P(16,1,4) = 1; P(11,6,4) = 1; P(15,5,4) = 1; P(14,9,4) = 1; 
P(4,4,4) = 1; P(7,7,4) = 1; P(10,10,4) = 1; P(13,13,4) = 1;
% others
P(:,:,5) = P(:,:,1)*P(:,:,4);
P(:,:,6) = P(:,:,1)*P(:,:,2);
P(:,:,7) = P(:,:,2)*P(:,:,4);


for k = 1:size(TableGBT,3)
    % merge if the same
    bases_ref = TableGBT(:,:,k);
    bases_ref = single(bases_ref);
    
    % compare if exactly the same
    cnt = 0;
    for col = 1:16
        if sum(bases(:,col) == bases_ref(:,col))==16 || sum(bases(:,col) == -bases_ref(:,col))==16
            cnt = cnt+1;
        end
    end
    if cnt == 16
        flag = 1;
        freq(k) = freq(k) + 1;
        return;
    end
    % merge if permutated
    partt = unique(bases(:,1));
    partt_ref = unique(bases_ref(:,1));
    % check: number of DC components & number of vertices within each
    % component
    if length(partt) == length(partt_ref) && sum(partt == partt_ref) == length(partt)
%         if length(partt) == 1
%             % different connected graphs
%             continue;
%         end
        % check mirror permutation
        for m = 1:7
            if sum(abs(bases(:,1)) == P(:,:,m)*abs(bases_ref(:,1)))==16 && sum(abs(bases(:,2)) == P(:,:,m)*abs(bases_ref(:,2))) == 16  % consider the first two DC bases
                flag = 1;
                flag_permute = 1;  % to think: transmit the permutation matrix or index?
                freq(k) = freq(k) + 1;
%                 if m > 4
%                     check = 1;
%                 end
                return;
            end
        end
    end  
    
end

% add
if flag == 0
    TableGBT(:,:,end+1) = bases;
    freq(end+1) = 1;
end

end