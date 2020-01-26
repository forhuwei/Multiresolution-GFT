function img2eps(file, name)
%------------------------------------
% Powered by leoman@ust.hk
% you can also specify the resolution
% print([filename(1:end-4),'.eps'],'-depsc','-r300');
%------------------------------------

if min(size(file)) == 1                            % input image name
   img = imread(file);
   imshow(img,'Border','tight',...                % display in a figure window without
          'InitialMagnification',100);            % a border at full magnification
   print([file(1:end-4),'.eps'],'-depsc');    % print the figure as eps

elseif nargin == 2 && ischar(name)
   img = file;                                    % input image
   imshow(img,'Border','tight',...                % display in a figure window without
           'InitialMagnification',100);           % a border at full magnification
   print([name,'.eps'],'-depsc');    % print the figure as a B&W eps
else
   disp('error parameter!');
end

disp('done!');
%close(gcf);
end