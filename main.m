clc
clear 

% number of agents
% nCk = [50, 70, 90, 100];
nCk = [100, 130, 170, 200];
% nCk = [20,20,20,20];
% nCk = [3, 4, 6, 7];
n = sum(nCk);
nx = 2;
n_1 = nCk(1);
n_2 = nCk(2);
n_3 = nCk(3);
n_4 = nCk(4);

%% Data load
% load("Data-16-05-2022/A_1.mat");
% load("Data-16-05-2022/A_2.mat");
% load("Data-16-05-2022/A_3.mat");
% load("Data-16-05-2022/A_4.mat");
load("Data-16-05-2022/init.mat");
load("Data-16-05-2022/Lext.mat");
load("Data-16-05-2022/Lint1.mat");
load("Data-16-05-2022/Lint2.mat");
load("Data-16-05-2022/Lint3.mat");
load("Data-16-05-2022/Lint4.mat");
    

%% State matrices

% Generates random 2d rotational matrices based on the number of agents in
% the network
A_1 = rotation_2d_matrix(nCk(1));
A_2 = rotation_2d_matrix(nCk(2));
A_3 = rotation_2d_matrix(nCk(3));
A_4 = rotation_2d_matrix(nCk(4));



%% Initial Value Generators
% init =  InitialConditions(nCk);


%% Generate random internal digraphs
% Generates random Laplacian of the weighted digraphs for clusters
% Lint1 = rand_lap(n_1);
% Lint2 = rand_lap(n_2);
% Lint3 = rand_lap(n_3);
% Lint4 = rand_lap(n_4);

% Generate the random directed internal digraphs for clusters
% Lint1 = weighted_directed_Lap_ran(n_1);
% Lint2 = weighted_directed_Lap_ran(n_2);
% Lint3 = weighted_directed_Lap_ran(n_3);
% Lint4 = weighted_directed_Lap_ran(n_4);

% Blockdiagonal internal Lapalcian of the network
Lint = blkdiag(Lint1, Lint2, Lint3, Lint4);

%% Parameters associated with the Laplacian matrix.

% This function returns the following
% PosLapEig: positive eigenvalues of the internal Laplacian matrices
% LambdaMin: minimum eigenvalue from each cluster
% ZeroIndex: index of the zero eigenvalue from each cluster. 
% LapNorm: Norm of the internal Laplacian Matrices.
[PosLapEig, LambdaMin, ZeroIndex, LapNorm] = LaplacianParameters(Lint, nCk);

%% External Laplacian Matrix
% Generates the random connected external graphs for the network.
% External Graph
Lext = laplacianER(nCk, 0.2,2);
% Lext = wt_dir_Lap (nCk);

% Network Laplacian
L = Lint + Lext;

%% Transformation Matrices
% Intra-cluster transformation
[H_l,H_r,Z_l,Z_r, H_l_bar,H_r_bar] = IntraClusterTranformationMatrices(nCk, Lint);

% Inter-cluster transformation
Ls = H_l_bar*Lext*H_r_bar;
[w_l, w_r, W_l, W_r, PosLapEigExt] = InterClusterTransformationMatrices(Ls);


%% Network Parameters
sigma_E = 200;  % External Gain
sigma_I = 500;  % Internal gain

% sigma_E = LambdaMin*sigma_I*(1/40) %Uncomment this line to make the parameter epsilon_1 = epsilon 2 and comment line with Sigma_E.   
% sigma_I = 50*sigma_E/min(Lint_norm);

gamma = sigma_I/sigma_E; % Ratio between the internal and the external gains. 

%% Network Parameters
norm_Lext = norm(Lext)
norm_Lint = norm(Lint)
mu = norm(Lext)/(gamma*LambdaMin)  % network parameter mu
epsilon_1 = 1/(sigma_E*norm(Ls))    % External perturbation parameter
epsilon_2 = epsilon_1*mu  % Internal perturnation parameter. 

%% Network Dynamics
% t_span = [0 1];
t_span = 0:0.00001:5;
A = blkdiag(A_1,A_2,A_3,A_4);
A_cl = (A- sigma_I*kron(Lint, eye(2)) - sigma_E*kron(Lext,eye(2)));
[t,x] = ode23(@(t,x) A_cl*x, t_span, init);
[x_Comp1, x_Comp2] =  StateComponents(x); % Separates the first and second components of the state x


%% Weighted Average Dynamics TSM1
A_s = H_l*blkdiag(A_1, A_2,A_3, A_4)*H_r;
A_agg = (A_s- sigma_E*H_l*kron(Lext,eye(nx))*H_r);

%Change initial conditions
x0_avg = H_l*init;
[t1,y]=  ode23s(@(t1,y) A_agg*y, t_span, x0_avg);
[y_Comp1, y_Comp2] =  StateComponents(y);

%% Intra-cluster error dynamics
xi_0 = Z_l*init;
P_33_3 = - (1/epsilon_2)*kron(diag(PosLapEig), eye(nx))/LambdaMin;
[t_xi, xi] = ode23s(@(t_xi,xi) P_33_3*xi, t_span, xi_0);
[xi_Comp1, xi_Comp2] =  StateComponents(xi);


%% TSM 2
% Emergent dynamics
P_11 = w_l*A_s*w_r;
x0_e = w_l*x0_avg;
[t_e, x_e] = ode23(@(t_e,x_e) P_11*x_e, t_span,x0_e);


% Inter-cluster synchronization dynamics
eta_0  = W_l*x0_avg;
P_22_2 = -  kron(PosLapEigExt,eye(nx))*(1/epsilon_1);
[t_eta, eta] = ode23s(@(t_eta,eta) P_22_2*eta, t_span, eta_0);
[eta_Comp1, eta_Comp2] =  StateComponents(eta);

%% Plots- Error 
figure 
IndexOfInterest_SemiLog = (t<=1) & (t>=0);
pl_intra_err = semilogx(t_xi(IndexOfInterest_SemiLog), xi_Comp1(IndexOfInterest_SemiLog,:));
hold on
pl_inter_err = semilogx(t_eta(IndexOfInterest_SemiLog),eta_Comp1(IndexOfInterest_SemiLog,:),'Color','#0072BD','LineStyle',':','LineWidth',1.5);
hold on 
pl_emerg = semilogx(t_e(IndexOfInterest_SemiLog),x_e(IndexOfInterest_SemiLog,1),'Color','#A2142F','LineStyle','--','LineWidth',1.5);
grid minor
legend ([pl_intra_err(1), pl_inter_err(1), pl_emerg], {'Intra-cluster Error Dyn (\xi_{f}(1))','Inter-cluster Error Dyn (\eta_{f}(1))','Emergent Dynamics (x_{s}(1))'})
xlabel('Time(s)')


%% Plots
% Long term state trajectories.
    indx_interest = (t <= 0.05) & (t >= 0);
indx_interest1 = (t_e <= 0.05) & (t_e >=0);
figure
plot(t(indx_interest), x_Comp1(indx_interest,:));
hold on 
p2 = plot(t_e(indx_interest1), x_e(indx_interest1,1), 'Color','#A2142F','LineStyle','--','LineWidth',2);
grid minor
xlabel('Time(s)')
legend ([p2], 'Emergent Dynamics (x_{s}(1))')
ylabel('x(1)')

% axes('position',[.58 .18 .28 .25])
% box on % put box around new pair of axes
% indexOfInterest = (t < 0.02) & (t >=0); % range of t near perturbation
% plot(t(indexOfInterest),x_Comp1(indexOfInterest,:)); % plot on new axes
% ylim([-200 250])
% axis fill
% grid minor

%% Plots
% Fast dynamics
% indx_interest = (t < 0.05) & (t >= 0);
% figure
% plot(t(indx_interest), x_Comp1(indx_interest,:));
% hold on 
% p2 = plot(t_e(indx_interest), x_e(indx_interest,1), 'Color','#A2142F','LineStyle','--','LineWidth',2);
% grid minor
% xlabel('Time(s)');
% legend ([p2], 'Emergent Dynamics (x_{e}(1))')
% ylabel('x(1)');
%     
% axes('position',[.60 .20 .28 .20])
% box on % put box around new pair of axes
% indexOfInterest = (t < 0.0003) & (t >=0); % range of t near perturbation
% plot(t(indexOfInterest),x_Comp1(indexOfInterest,:)); % plot on new axes
% ylim([-200 250])
% axis fill
% grid minor

%% Plots
% Plots log vs x1
indx_interest = (t < 1) & (t >= 0);
figure
p1 = semilogx(t(indx_interest), x_Comp1(indx_interest,:));
hold on 
p2 = semilogx(t_e(indx_interest), x_e(indx_interest,1), 'Color','#A2142F','LineStyle','--','LineWidth',2);
grid minor
xlabel('Time (s)')
ylabel('x(1)')
legend ([p2], 'Emergent Dynamics (x_{s}(1))')

%% Plots 
% x1 vs x2
% figure
% l = 20;
% p1 = plot(x_1(l:end,:), x_2(l:end,:));
% grid minor
% xlabel('x(1)')
% ylabel('x(2)')







