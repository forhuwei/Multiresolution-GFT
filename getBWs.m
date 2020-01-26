function [ BW_s ] = getBWs( BWb,hX,wX,L )

BW_s = zeros(2*hX,2*wX);
for i = 1:2*hX
    for j = 1:2*wX
        if i==38 && j==326
            flag=1;
        end
       if mod(i,2) == 0 && mod(j,2) == 0  % diagonal edge
           ix = L*(i-2)+2;
           iy = L*(j-2)+2;
           for k = 1:L 
               if BWb(ix,iy) == 1
                   BW_s(i,j) = 1;
                   break;
               end
               ix = ix+2; 
               iy = iy+2;
           end
           if BW_s(i,j) == 0
               ix = L*i;
               iy = L*(j-2)+2;
               for k = 1:L 
                   if BWb(ix,iy) == 1
                       BW_s(i,j) = 1;
                       break;
                   end
                   ix = ix-2; 
                   iy = iy+2;
               end
           end
       end
       if mod(i,2) == 0 && mod(j,2) ~= 0  % vertical edge
           ix = L*(i-2)+2;
           iy = L*(j-1)+1;
           for k = 1:L 
               if BWb(ix,iy) == 1
                   BW_s(i,j) = 1;
                   break;
               end
               ix = ix+2;    
           end
       end
       if mod(i,2) ~= 0 && mod(j,2) == 0  % horizontal edge
           iy = L*(j-2)+2;
           ix = L*(i-1)+1;
           for k = 1:L
               if BWb(ix,iy) == 1
                   BW_s(i,j) = 1;
                   break;
               end
               iy = iy+2;    
           end
       end
    end
end


end

