function L = wt_dir_Lap (nCk)
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



% nCk = [3, 4, 6, 7];
% n=  5;
% density = 1;

% density = 1e+3;
n = sum(nCk);
Adj = 5*sprand(n,n,0.0005);
Adj = full(Adj);
Adj = abs(Adj);

Adj_0 = ones(n,n);

for i = 1:length(nCk)
   if i == 1
        st = 1;
        ed = nCk(1);
        for j = st:ed
            Adj_0(st:ed,j)=0;
        end
   else
       st = sum(nCk(1:i-1))+1;
       ed = st + nCk(i)-1;
       for j = st:ed
            Adj_0(st:ed,j)=0;
       end
   end


end

Adj = floor(Adj.*Adj_0);

Deg = diag(Adj*ones(n,1));
L = Deg- Adj;
end