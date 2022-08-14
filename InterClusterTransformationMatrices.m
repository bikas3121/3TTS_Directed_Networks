function [w_l, w_r, W_l, W_r, D] = InterClusterTransformationMatrices(Ls)


% Ls = rand_lap(5);

% Calculate the eigenvalues and the eigenvector
[V,D] = eig(Ls);
InV = inv(V);

% Select the real parts only
V = real(V);
InV = real(InV);
D = real(D);
D_v = diag(D);

% Index of zero eigenvalue
for i = 1:length(D_v)
    if D_v(i) <= 0.000001
        index_zero = i;
    end
end

% Eigenvalue matrix without zero eigenvalue
D_v(index_zero) = [];
D = diag(D_v);

% Scaling of the eigenvectors
fac = V(index_zero,index_zero);
V = V*(1/fac);
InV = InV*fac;

% Tranformation vectors
% Left eigenvectors
w_l = kron(InV(index_zero,:),eye(2));
W_l = InV;
W_l(index_zero,:) = [];
W_l = kron(W_l,eye(2));
% Right Eigenvectors
w_r = kron(V(:,index_zero),eye(2));
W_r = V;
W_r(:,index_zero) = [];
W_r = kron(W_r,eye(2));

end