function [bCount] = arithenco_mine(seq,prob)

b = 0;   % interval base
l = 1;   % interval length
bCount = 0;  
gamma = 2;  % interval rescaling factor

for ind = 1:length(seq),
    
    dcode = seq(ind);
    
    % update interval 
    b = b + sum(prob(1:dcode-1)) * l;
    l = prob(dcode) * l;
    
    % carry propagation
    if b>=1  % there is a carry
        b = b-1;
        d = Propagate_Carry(bCount,d);
    end
    
    % interval renormalization
    while( l <= 0.5 )
        bCount = bCount+1;
        l = 2*l;
        if b >= 0.5
            d(bCount) = 1;  % output bit 1
            delta = 0.5;
            b = gamma*(b-delta);
        else
            d(bCount) = 0;  % output bit 0
            delta = 0;          
            b = gamma*(b-delta);
        end
    end
end

% step 3: choose final code value
bCount = bCount+1;
if b <= 0.5
    d(bCount) = 1;
else
    d(bCount) = 0;
    d = Propagate_Carry(bCount-1,d);
end

% remove the 0's at the end
% while(d(end)==0)
%     d(end) = [];
%     bCount = bCount-1;
% end

end
