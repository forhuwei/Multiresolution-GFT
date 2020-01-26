function [ d ] = Propagate_Carry( t,d )

n = t;

% complement all the outstanding bits until the first 0-bit is complemented
while (d(n) == 1)  
    d(n) = 0;  
    n = n - 1;
    if n==0
        break;
    end
end
if n>=1
    d(n) = 1; 
else
    len = length(d);
    x=zeros(len+1,1);
    x(1) = 1;
    x(2:end) = d;
    d=x;
end

end

