% *************************************************************
% author:   Wei Hu                                           **
% date:     09/01/2012                                       **
% modified: 12/20/2012                                       **
% purpose:  Construct a s,t-graph for Hochbaum's algorithm   **
% *************************************************************

function [ value,cut ] = MRF_hpf( N1,N,bSize,dep,rho,alpha )

    f = dep(:);
    
    % construct s,t-graph 
    P = zeros(N+2); 

    % calculate the subgradient of Gi() at alpha = 5
    G_i_d = zeros(N,1);
    for i = 1:N1
        G_i_d(i) = 0.13*(f(i)-f(i+bSize))^2 - 0.13*alpha*rho;
    end

    k = 1;
    for i = (N1+1):N
        G_i_d(i) = 0.13*(f(k)-f(k+1))^2 - 0.13*alpha*rho;
        if mod(k,bSize) == bSize-1 
            k = k+2;
        else
            k = k+1;
        end
    end

    % arc to source: max(0, G_v_d(alpha))
    for j = 2:(N+1)
        P(1,j) = max(0,G_i_d(j-1));
        P(j,1) = P(1,j);
    end

    % arc to sink: max(0, -G_v_d(alpha))
    for j = 2:(N+1)
        P(N+2,j) = max(0,-G_i_d(j-1));
        P(j,N+2) = P(N+2,j);
    end

    % arc among nodes: u_ij
    % horizontal graph edges 
    neighbor = zeros(24,6);
    neighbor(1,1:3) = [2 13 16];
    neighbor(5,1:3) = [6 16 19];
    neighbor(9,1:3) = [10 19 22];
    neighbor(4,1:3) = [3 15 18];
    neighbor(8,1:3) = [7 18 21];
    neighbor(12,1:3) = [11 21 24];
    
    neighbor(2,:) = [1 3 13 14 16 17];
    neighbor(6,:) = [5 7 16 17 19 20];
    neighbor(10,:) = [9 11 19 20 22 23];
    neighbor(3,:) = [2 4 14 15 17 18];
    neighbor(7,:) = [6 8 17 18 20 21];
    neighbor(11,:) = [10 12 20 21 23 24];
    
    neighbor(13,1:3) = [1 2 16];
    neighbor(14,1:3) = [2 3 17];
    neighbor(15,1:3) = [3 4 18];
    neighbor(22,1:3) = [9 10 19];
    neighbor(23,1:3) = [10 11 20];
    neighbor(24,1:3) = [11 12 21];
    
    neighbor(16,:) = [13 19 1 2 5 6];
    neighbor(17,:) = [14 20 2 3 6 7];
    neighbor(18,:) = [15 21 3 4 7 8];
    neighbor(19,:) = [16 22 5 6 9 10];
    neighbor(20,:) = [17 23 6 7 10 11];
    neighbor(21,:) = [18 24 7 8 11 12];
        
    for i = 1:N
        for j = 1:6
            if neighbor(i,j) ~= 0
                P(i+1,neighbor(i,j)+1) = rho;
            end
        end
    end

    % min-cut
    capacity = sparse(P);            
    source = 1;
    sink = N+2;
    [value,cut] = hpf(capacity,source,sink);
    
%     if cut(2:end) == 0  % all weights 1
%         return;
%     end
%     for k = 1:4
%         % update V
%         ind = find(cut(2:end-1)==1);
%         for i = 1:length(ind)
%             P(ind(i)+1,:) = [];
%             P(:,ind(i)+1) = [];
%             ind = ind-1;
%         end
%         capacity = sparse(P); 
%         sink = sink - length(ind);
%         [value,cut] = hpf(capacity,source,sink);
%         if cut(2:end) == 0  % all weights 1
%             return;
%         else
%             check = 1;
%         end
%     end
end

