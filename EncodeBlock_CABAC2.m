function [rate, rate_edges] = EncodeBlock_CABAC2(Yb,BWb,QP)

Yb = Yb(:);
rate = 0;
rate_edges = 0;

fin = 'coeff_input_block.txt';
fin2 = 'coeff_input_edges.txt';
fout = 'coeff_output_block.txt';
frate = 'rate_output_block.txt';
frate2 = 'rate_output_edges.txt';

fid = fopen(fin2,'w');
for i = 1:length(BWb(:))-1
    fprintf(fid, '%d, ', BWb(i));
end
fprintf(fid, '%d', BWb(length(BWb)));
fclose(fid);

fid = fopen(fin,'w');
fprintf(fid, '%d, ', length(Yb(:)));
for i = 1:length(Yb(:))
    fprintf(fid, '%d, ', Yb(i));
end
fprintf(fid, '%d', Yb(length(Yb(:))));
fclose(fid);

temp = 'lencod_oneblock2.exe -f encoder.cfg > results.txt';
system(temp);

% Extract bitrate
fid = fopen(frate,'r');
rate = fscanf(fid, '%d', 1);
fclose(fid);

% Extract bitrate
fid = fopen(frate2,'r');
rate_edges = fscanf(fid, '%d', 1);
fclose(fid);