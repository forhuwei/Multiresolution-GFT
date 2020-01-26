% *************************************************************
% author:   Gene Cheung                                      **
% date:     10/28/2012                                       **
% modified: 01/08/2013                                       **
% purpose:  arithmetic edge encoding                         **
% *************************************************************

% given direction code (1-4 means North, East, South, West),
% compute AEC codeword (d).

% reference: Amir Said, "Introduction to Arithmetic Coding - 
% Theory and Practice", April 2004. 
% http://www.hpl.hp.com/techreports/2004/HPL-2004-76.pdf

function [d, bCount] = aec (dcode)

% step 0: initialize var
fitLength = 3;    % length of input vector for bestFitLine
b = 0;   % interval base
l = 1;   % interval length
bCount = 0;  
gamma = 2;  % interval rescaling factor

for ind = 1:length(dcode),
    
    % step 1: find best estimate of prob of 3 poss. directions
    % (must be same at decoder in aed.m)

    if ind == 1,
        prob = ones(1,4) / 4.0;    % no end-code for first chain
    else
        if ind > fitLength,
%             rho = 4;
            prob = bestFitLine(dcode(ind-fitLength:ind-1));
        else
%             rho = 2;
%             prob = bestFitLine(dcode(1:ind-1));
            prob = ones(1,3) / 3.0;
        end
        % add end-code prob, can be any function of ind
        endProb = 0;   % here end-code prob unused, code the symbol number instead
        prob = [prob*(1.0-endProb)  endProb];
    end
    % ****** (must be same at decoder in aed.m) ******

    % step 2: actual arithmetic coding:  given prob. vector, index,
    % compute new interval

    if ind == 1,
        ccode = dcode(1);  % fist chain is in absolute term  
    elseif ind <= length(dcode),
        if dcode(ind) == dcode(ind-1),
            ccode = 2;     % straight
        elseif dcode(ind) - dcode(ind-1) == 1,
            ccode = 3;     % right
        elseif dcode(ind) ==1 && dcode(ind-1) == 4,
            ccode = 3;
        else
            ccode = 1;     % left
        end
    else
        ccode = 4;         % end-code
    end
    
    % update interval 
    b = b + sum(prob(1:ccode-1)) * l;
    l = prob(ccode) * l;
    
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
while(d(end)==0)
    d(end) = [];
    bCount = bCount-1;
end


end

