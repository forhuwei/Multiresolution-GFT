function [ psnr ] = calPSNR( I, I_target )

I = double(I);
I_target = double(I_target);

[width height] = size(I);
% R1 = I(:,:,1);
% G1 = I(:,:,2);
% B1 = I(:,:,3);
% R2 = I_target(:,:,1);
% G2 = I_target(:,:,2);
% B2 = I_target(:,:,3);

mse = sum(sum(((I-I_target).^2)))/(width*height);
psnr = 10*log10(255*255/mse);

end

