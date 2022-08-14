function [H_l,H_r,Z_l,Z_r, H_l_bar,H_r_bar] = IntraClusterTranformationMatrices(nCk, Lint)

% nCk = [5, 4, 6, 7];
% 
% Lint1 = weighted_directed_Lap_ran(nCk(1));
% Lint2 = weighted_directed_Lap_ran(nCk(2));
% Lint3 = weighted_directed_Lap_ran(nCk(3));
% Lint4 = weighted_directed_Lap_ran(nCk(4));
% 
% Lint = blkdiag(Lint1, Lint2, Lint3, Lint4);

% Transformation Matrices

% Block diagonal transformation matrices with the kronecker product of the
% dimension of the system
H_l = [];
H_r = [];
Z_l = [];
Z_r = [];

% Block diagonal transformation matrices without the kronecker product of
% the dimension of the systems. 
H_l_bar = [];
H_r_bar = [];

% Matrix Generations
for i = 1:length(nCk)
    if i == 1      
        Lint_i = Lint(1:nCk(i), 1:nCk(i));

        [V_i,D_i] = eig(Lint_i); % Calculate the eigenvalue and right eigenvector of the Laplacian
        InV_i = inv(V_i); % Calculate the left eigenvectors of the Laplacian

        % Takes only the real parts of the D_i, V_i and InV_i
        V_i = real(V_i); 
        InV_i = real(InV_i);
        D_i = real(diag(D_i));
        
        % Determine the index of the eigenvalue that is zero
        for j = 1:length(D_i)
            if D_i(j) <= 0.0000001
                index_zero = j;
            end
        end
        
        % Scale the eigenvector  to make the eigenvectors
        % associated with the zero eigenvalue as vector of ones.
        fac_i = V_i(index_zero,index_zero);
         V_i = V_i*(1/fac_i);
         InV_i = InV_i*fac_i;

        % Left and Right Eigenvectors Associated with zero eigenvalues
         H_l_i = InV_i(index_zero,:);
         H_r_i = V_i(:,index_zero);

        % Left and right eigenvectors associated with positive eigenvalues
         Z_l_i = InV_i;
         Z_l_i(index_zero,:) =[];
         Z_r_i = V_i;
         Z_r_i(:,index_zero) =[];

         H_l = blkdiag(H_l, kron(H_l_i,eye(2)));
         H_r = blkdiag(H_r, kron(H_r_i,eye(2)));
         Z_l = blkdiag(Z_l, kron(Z_l_i,eye(2)));
         Z_r = blkdiag(Z_r, kron(Z_r_i,eye(2)));

         H_l_bar = blkdiag(H_l_bar, H_l_i);
         H_r_bar = blkdiag(H_r_bar, H_r_i);
    else
        
        st = sum(nCk(1:i-1))+1;
        ed = st+nCk(i)-1;
        Lint_i = Lint(st:ed, st:ed);

        [V_i,D_i] = eig(Lint_i); % Calculate the eigenvalue and right eigenvector of the Laplacian
        InV_i = inv(V_i); % Calculate the left eigenvectors of the Laplacian

        % Takes only the real parts of the D_i, V_i and InV_i
        V_i = real(V_i); 
        InV_i = real(InV_i);
        D_i = real(diag(D_i));
        
        % Determine the index of the eigenvalue that is zero
        for j = 1:length(D_i)
            if D_i(j) <= 0.0000001
                index_zero = j;
            end
        end
        
        % Scale the eigenvector  to make the eigenvectors
        % associated with the zero eigenvalue as vector of ones.
        fac_i = V_i(index_zero,index_zero);
         V_i = V_i*(1/fac_i);
         InV_i = InV_i*fac_i;

        % Left and Right Eigenvectors Associated with zero eigenvalues
         H_l_i = InV_i(index_zero,:);
         H_r_i = V_i(:,index_zero);

        % Left and right eigenvectors associated with positive eigenvalues
         Z_l_i = InV_i;
         Z_l_i(index_zero,:) =[];
         Z_r_i = V_i;
         Z_r_i(:,index_zero) =[];

         H_l = blkdiag(H_l, kron(H_l_i,eye(2)));
         H_r = blkdiag(H_r, kron(H_r_i,eye(2)));
         Z_l = blkdiag(Z_l, kron(Z_l_i,eye(2)));
         Z_r = blkdiag(Z_r, kron(Z_r_i,eye(2)));

         H_l_bar = blkdiag(H_l_bar, H_l_i);
         H_r_bar = blkdiag(H_r_bar, H_r_i);
    end
    
end
end