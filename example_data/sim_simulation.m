function data = sim_simulation(seed)
% Simulation for the paper with R = 10 and L = 1
%Output : Data generated from the models with parameters.
%Input : 
% seed : number of the random seed from which data was generated. Saved for
% replication purposes.
% L : number of lags
% G : number of groups
% n : number of subjects
% R : number of nodes
% T : number of time points
% Set seed:
rng(seed)

L = 1;
G = 2;
n = 20;
R = 10;
T = 400;
% Create group labels
eta = [ones(1,10),2*ones(1,10)];
% Define underlying structural connectivity
DTI_mat = [];
g_1 = [0.1  0   0   0   0    0   0   0   0   0 ;
       0.3 0.5  0   0   0    0   0   0   0   0 ;
       0.1 0.1 0.35 0   0    0   0   0   0   0 ;
       0.6 0.2 0.65 0.3 0    0   0   0   0   0 ;
       0.2 0.1 0.2 0.25 0.1  0   0   0   0   0 ;
       0.1 0.6 0.7 0.1  0.85 0.5 0   0   0   0 ;
       0.7 0.3 0.1 .2  0.3  0.3 0.1   0   0   0;
       0.3 0.1 0.5 .15  0.25  0.1 0.05   0.25   0   0;
       0.1 0.5 0.2 .4  0.1  0.25 0.12    0.3   0.1   0;
       0.4 0.15 0.4 .3  0.15  0.2 0.3   0.15   .3   0.2];
DTI_mat.one = g_1 + tril(g_1,-1).';
g_2 = [0.2  0   0   0   0    0   0   0   0   0 ;
       0.5 0.3  0   0   0    0   0   0   0   0 ;
       0.3 0.1 0.55 0   0    0   0   0   0   0 ;
       0.1 0.3 0.15 0.1 0    0   0   0   0   0 ;
       0.4 0.1 0.5 0.35 0.36  0   0   0   0   0 ;
       0.5 0.1 0.5 0.3  0.1 0.1 0   0   0   0 ;
       0.1 0.1 0.1 .6  0.2  0.1 0.7   0   0   0 ;
       0.2 0.15 0.5 .1  0.1  0.2 0.25   0.25   0   0 ;
       0.3 0.05 0.1 .4  0.25  0.5 0.4   0.3   0.1   0 ;
       0.1 0.3  0.4 .5  0.05  0.15 0.2   0.15   0.25   0.3 ];
       
DTI_mat.two = g_2 + tril(g_2,-1).';
% Probit regression sparsity parameter
alpha0_sim =[-1.5,-1.5]; % baseline
% Probit regression beta parameter %
alpha1_sim = [5,5]; % beta
% Sampling of non-zero effective connectivities using standard normal CDF
% probabilities using auxiliary variable z
z_sim = [];
DTI_vec = [];
gamma_sim =[];
names = {'one','two'};
for i = 1:G
DTI_vec.(names{i})=reshape(DTI_mat.(names{i}),R^2*L,1);
z_sim.(names{i}) =normrnd(alpha0_sim(i)+alpha1_sim(i)*DTI_vec.(names{i}),1,R^2,1);
% z_sim.(names{i}) =mvnrnd((alpha0_sim(i)+alpha1_sim(i)*DTI_vec.(names{i}))',eye(L*R^2),1);
gamma_sim.(names{i})=(z_sim.(names{i})>0)*1;
end
% Simulate group effective connectivity from unif(0,0.5) if non-zero
%gamma_mat =[];
PHI_g=[];
for i = 1:G
    PHI_g.(names{i}) = reshape(gamma_sim.(names{i}).*round(unifrnd(0,0.4,R^2,1),2),R,R);  
end
% Covariance of VAR error term
XI_sim = eye(R);

% Simulate subject level effective connectivities by adding subject-specific
% random matrix to group matrix, with eigenvalues (-0.4, -0.25, -0.10, 0.05, 
% 0.20)
% Create a string to store subject names
subjects = cell(n,1);
for i = 1:n
    subjects{i}=strcat('s_',num2str(i));
end
ev = [-0.4,-0.25,-0.1,0.05,0.2,-0.3, 0.1, 0.1,-0.3,-0.15];
phi_s=[];
beta_s =[];
beta_true_cat=[];
for i = 1:n
    [Q,R_comp]= qr(normrnd(0,1,R,R));
    ph =diag(R_comp)./abs(diag(R_comp));
    Z =transpose(Q*diag(ph))*diag(ev);
    phi_s.(subjects{i})=PHI_g.(names{eta(i)})+Z;
    beta_s.(subjects{i})=reshape(PHI_g.(names{eta(i)})+Z,R^2,1);
    beta_true_cat=[beta_true_cat; beta_s.(subjects{i})];
end
%Simulate subject-level time-series from VAR model of order 1
Y = zeros(T+1,R,n);
X = zeros(T,R,n);
Y(1,:,:)=0;
for s = 1:n
    for t = 2:T+1
        Y(t,:,s)=Y(t-1,:,s)*phi_s.(subjects{s})+mvnrnd(zeros(R,1),XI_sim);
        X(t-1,:,s)=Y(t,:,s);
    end
end
data =struct('Y',Y,'X',X,'DTI_vec',DTI_vec,'eta',eta,'gamma_sim',gamma_sim,...
    'beta_true_cat',beta_true_cat,'PHI_g',PHI_g);
end

