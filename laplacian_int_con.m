function [Lint] =  laplacian_int_con(n,p,Kreg)
% 
% p =0.1;
% n = 5;
% Kreg  = 4;
[G]=erdosRenyi(n,p,Kreg);

adj = full(G.Adj);

% make adj symmetric
for  i =1:n
    for j = 1:n
        if adj(i,j) == 1
            adj(j,i) = 1;
        end
    end
    adj(i,i) =0;
end


deg = diag(adj*ones(n,1));
Lint = deg-adj;
end