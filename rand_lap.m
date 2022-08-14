function L = rand_lap (n)
% This functio rand_lap generates Laplacian for the weighted directed
% graph. The weihts vary between the 1-100 that can be changed as per the
% requirement. 


 

rn = randi(2,n,n);
Adj = rn; % define the random matrix as weighted adjacency matrix

% remove the zeros from the diagonal
for i = 1:n
    Adj(i,i) =0;
end
% Adj

% generate the indegree matrix. 
Deg = zeros(n,n);
for i = 1:n
    Deg(i,i) = sum(Adj(i,:));
end
% Deg

L = Deg-Adj;

end