function [bCount] = AriCod(coeff)

unix = unique(coeff);
if length(unix) == 1
    bCount = 0;
    return;
end

p = zeros(size(unix));
for i = 1:length(coeff)
    for k = 1:length(unix)
        if coeff(i) == unix(k)
            p(k) = p(k) + 1;
            seq(i) = k;
            continue;
        end
    end
end

p = p/sum(p);
bCount = arithenco_mine(seq,p);

end