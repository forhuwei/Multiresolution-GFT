function [ diff ] = computeDiff( f )

diff = zeros(1,24);
for i = 1:12
    diff(i) = (f(i)-f(i+4))^2;
end

k = 1;
for i = 13:24
    diff(i) = (f(k)-f(k+1))^2;
    if mod(k,4) == 4-1 
        k = k+2;
    else
        k = k+1;
    end
end

end

