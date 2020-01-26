function [ Rate ] = proxyTrRate( EdgeW )
% rate proxy of weak edges

Rate = 1;

for j = 2:2:6
    for i=1:2:5
        Rate = Rate + abs(EdgeW(i,j)-EdgeW(i+2,j));
    end
end

for j = 1:2:5
    for i=2:2:6
        Rate = Rate + abs(EdgeW(i,j)-EdgeW(i,j+2));
    end
end


end

