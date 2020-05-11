function data = sim_simulation_lag_2(seed)
rng(seed)
L = 2;
G = 2;
n = 20;
R = 10;
T = 400;
% Create group labels
eta = [ones(1,10),2*ones(1,10)];
% Define underlying structural connectivity
DTI_mat = [];
% Create Structrual Prior
vals = unifrnd(0.3,0.7,R/2*(R + 1),1);
vals(randsample(1:R/2*(R + 1),35)) = 0.1;
B = zeros(10); % define B of required size
[ii jj] = ndgrid(1:10); % ii and jj are row and column indices respectively
B(ii>=jj) = vals; 
mat_1 = B + tril(B,-1).';
vals = unifrnd(0.3,0.7,R/2*(R + 1),1);
vals(randsample(1:R/2*(R + 1),35)) = 0.1;
B = zeros(10); % define B of required size
[ii jj] = ndgrid(1:10); % ii and jj are row and column indices respectively
B(ii>=jj) = vals; 
mat_2 = B + tril(B,-1).';
DTI_mat.one  = [mat_1;mat_2];
vals = unifrnd(0.3,0.7,R/2*(R + 1),1);
vals(randsample(1:R/2*(R + 1),35)) = 0.1;
B = zeros(10); % define B of required size
[ii jj] = ndgrid(1:10); % ii and jj are row and column indices respectively
B(ii>=jj) = vals; 
mat_1 = B + tril(B,-1).';
vals = unifrnd(0.3,0.7,R/2*(R + 1),1);
vals(randsample(1:R/2*(R + 1),35)) = 0.1;
B = zeros(10); % define B of required size
[ii jj] = ndgrid(1:10); % ii and jj are row and column indices respectively
B(ii>=jj) = vals; 
mat_2 = B + tril(B,-1).';
DTI_mat.two  = [mat_1;mat_2];
% Make diagonals non-zero
for i = 1:10
    for j = 1:2
        DTI_mat.one(i+(j-1)*10,i) = min(DTI_mat.one(i+(j-1)*10,i)+0.5,1);
        DTI_mat.two(i+(j-1)*10,i) = min(DTI_mat.two(i+(j-1)*10,i)+0.5,1);
    end
end

% Probit regression sparsity parameter
alpha0_sim =[-1.5,-1.5]; % baseline
% Probit regression beta parameter %
alpha1_sim = [3,3]; % beta
% Sampling of non-zero effective connectivities using standard normal CDF
% probabilities using auxiliary varaible z
z_sim = [];
DTI_vec = [];
gamma_sim =[];
names = {'one','two'};
for i = 1:G
DTI_vec.(names{i})=reshape(DTI_mat.(names{i}),R^2*L,1);
z_sim.(names{i}) =normrnd(alpha0_sim(i)+alpha1_sim(i)*DTI_vec.(names{i}),1,R^2*L,1);
gamma_sim.(names{i})=(z_sim.(names{i})>0)*1;
end

PHI_g=[];
for i = 1:G
    PHI_g.(names{i}) = reshape(gamma_sim.(names{i}).*round(unifrnd(-0.4,0.4,R^2*L,1),2),R*L,R);  
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

ev = unifrnd(-0.2,0.1,1,10);
phi_s=[];
beta_s =[];
beta_true_cat=[];
for i = 1:n
    Z = [];
    for j = 1:2
        [Q,R_comp]= qr(normrnd(0,1,R,R));
        ph =diag(R_comp)./abs(diag(R_comp));
        temp =transpose(Q*diag(ph))*diag(ev);
        Z = [Z;temp];
    end
    phi_s.(subjects{i})=PHI_g.(names{eta(i)})+Z;
    beta_s.(subjects{i})=reshape(PHI_g.(names{eta(i)})+Z,R^2*L,1);
    beta_true_cat=[beta_true_cat; beta_s.(subjects{i})];
end
%Simulate subject-level time-series from VAR model of order 1
Y = zeros(T+2,R,n);
X = zeros(T,R,n);
Y(1,:,:)=0;
Y(2,:,:)=0.1;
for s = 1:n
    for t = 3:T+2
        Y(t,:,s) = [Y([t-1],:,s),Y([t-2],:,s)]*phi_s.(subjects{s})+mvnrnd(zeros(R,1),XI_sim);
        X(t-2,:,s)=Y(t,:,s);
    end
end
data =struct('Y',Y,'X',X,'DTI_vec',DTI_vec,'eta',eta,'gamma_sim',gamma_sim,...
    'beta_true_cat',beta_true_cat,'PHI_g',PHI_g);

end