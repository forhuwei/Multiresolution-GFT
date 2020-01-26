function [ err ] = LookupErr( BW,bases )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

w = ones(24,1);
[indX,indY] = find(BW==1);
for k=1:length(indX)
    if mod(indX(k),2) ~= 0
        nodeX = (indX(k)+1)/2;
        nodeY = indY(k)/2;
        node1 = 4*(nodeY-1) + nodeX;
        w(node1) = 0;
    else
        nodeX = indX(k)/2;
        nodeY = (indY(k)+1)/2;  
        node1 = 4*(nodeY-1) + nodeX;
        if node1<=3
            w(node1+12) = 0;
        elseif node1 <= 7
            w(node1+11) = 0;
        elseif node1 <= 11
            w(node1+10) = 0;
        else
            w(node1+9) = 0;
        end
            
    end
end
    
basis = bases(:,1);
diff = zeros(24,1);
for k=1:12
    diff(k) = basis(k)-basis(k+4);
end
i=1;
for k=13:24
    diff(k) = basis(i)-basis(i+1);
    if mod(i,4) == 3
        i = i+2;
    else
        i = i+1;
    end
end

basis2 = bases(:,2);
diff2 = zeros(24,1);
for k=1:12
    diff2(k) = basis2(k)-basis2(k+4);
end
i=1;
for k=13:24
    diff2(k) = basis2(i)-basis2(i+1);
    if mod(i,4) == 3
        i = i+2;
    else
        i = i+1;
    end
end

basis3 = bases(:,3);
diff3 = zeros(24,1);
for k=1:12
    diff3(k) = basis3(k)-basis3(k+4);
end
i=1;
for k=13:24
    diff3(k) = basis3(i)-basis3(i+1);
    if mod(i,4) == 3
        i = i+2;
    else
        i = i+1;
    end
end

    
err1 = sum(w.*(diff.^2));
err2 = sum(w.*(diff2.^2));
err3 = sum(w.*(diff3.^2));
% err = max([err1,err2,err3]);
err = err1+err2+err3;

end

