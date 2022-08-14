function [A] = rotation_2d_matrix(nk)
% The function rotation_2d_matrix generates the 2-dimensional rotational
% matrix and put them in block diagonal form with nk representing the
% number of blocks. 

% Input:
%nk - number of blocks 

% Output:
% A - block diagonal matrix with dimension 2nk times 2nk


fac = 6000*rand([1,nk]);
A = [];
for i = 1:length(fac)
    A_k = fac(i)*[0,1;-1,0];
    A = blkdiag(A, A_k);
end

end