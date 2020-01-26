% *************************************************************
% author:   Wei Hu                                           **
% date:     12/20/2012                                       **
% modified: 12/20/2012                                       **
% purpose:  construct the adjacency matrix for WGBT          **
% *************************************************************

function [W,FlagW, EdgeW] = AdjMatW(cut,bSize,weight,img,T)

FlagW = 0;
EdgeW = zeros(bSize*2);
W =     [0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0;
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

ind = find(cut(2:end) == 1);
for k = 1:length(ind)
    ind_k = ind(k);
    if ind_k <= 12
        node1 = ind_k;
        node2 = node1+4;
    elseif ind_k <= 15
        node1 = ind_k-12;
        node2 = node1+1;
    elseif ind_k <= 18
        node1 = ind_k-11;
        node2 = node1+1;
    elseif ind_k <= 21
        node1 = ind_k-10;
        node2 = node1+1;
    else
        node1 = ind_k-9;
        node2 = node1+1;
    end
    
    diff = abs(img(node1)-img(node2));
    
    if diff <= T
        W(node1,node2) = weight;
        W(node2,node1) = weight;
        if node2 == node1+4
            i = mod(node1,bSize);
            if i == 0
                i = bSize;
            end
            j = (node1-i)/bSize+1;
            i = 2*i-1;
            j = j*2;
        else
            i = mod(node1,bSize); 
            if i == 0
                i = bSize;
            end
            j = (node1-i)/bSize+1;
            i = 2*i;
            j = 2*j-1;
        end
        EdgeW(i,j) = 1;
        FlagW = 1;  
    end
    
    if diff > T
        W(node1,node2) = 0;
        W(node2,node1) = 0;
        if node2 == node1+4
            i = mod(node1,bSize);
            if i == 0
                i = bSize;
            end
            j = (node1-i)/bSize+1;
            i = 2*i-1;
            j = j*2;
        else
            i = mod(node1,bSize); 
            if i == 0
                i = bSize;
            end
            j = (node1-i)/bSize+1;
            i = 2*i;
            j = 2*j-1;
        end
        EdgeW(i,j) = 1;
        FlagW = 0;  
    end
end

end





        
