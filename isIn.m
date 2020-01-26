function [ flag ] = isIn( array1,array2 )
% check if array1 (single row) is in array2 or not
% return 1 if is in 

if isempty(array1) || isempty(array2)
    flag = 0;
    return;
end

[rows cols] = size(array2);
flag = 0;

for i = 1:rows
    if sum(array1 == array2(i,:)) == cols
        flag = i;
        return;
    end
end

end

