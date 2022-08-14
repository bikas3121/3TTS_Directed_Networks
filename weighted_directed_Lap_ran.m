function [L] = weighted_directed_Lap_ran (n)
% The function wt_dir_Lap generates the external laplacian of the weighted directed
% graph. The inputs are;
%   n = number of nodes in the graph. 
%   density = density of the edges in the graphs, roughly. Value should be
%   between [0,1].
% The output is the external weighted directed Laplacian matrix. The
% external Laplacian does not have the connections inside the clusters. A
% sparse adjacency matrix is randomly generated and the connections inside the
% cluster are removed. Then the Laplacian generated from the corresponding
% adjacency matrix. 




% n=  6;
% density = 1;

% density = 1e+3;
% n = sum(nCk);
Adj = 7*sprand(n,n,0.3);
Adj = full(Adj);
Adj = abs(Adj);
Adj = floor(Adj);
% rn = randi(5,n,n);

% Remove diagonal elements (i.e. no self-loops)
for i = 1:n
    Adj(i,i) =0;
end
Deg= diag(Adj*ones(n,1));
L = Deg-Adj;


%% Alternative using the inbuilt function 'randi'

% %r = randi([a,b],m,n) 
% % The function 'randi' creates the random matrix of size (m,n) with values
% % between  (a,b).
% 
% rn = randi([0,3],n,n);
% Adj = rn; % define the random matrix as weighted adjacency matrix
% 
% % remove the zeros from the diagonal
% for i = 1:n
%     Adj(i,i) =0;
% end
% 
% Deg= diag(Adj*ones(n,1));
% 
% 
% % Deg
% 
% L = Deg-Adj;

end