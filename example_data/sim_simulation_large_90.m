function data = sim_simulation_large_90(seed)
rng(seed)
L = 1;
G = 2;
n = 100;
%n = 200;
R = 90;
T = 150;
% Create group labels
eta = [ones(1,50),2*ones(1,50)];
%eta = [ones(1,100),2*ones(1,100)];
% Define underlying structural connectivity
DTI_mat = [];
% Create Structrual Prior
vals = unifrnd(0.3,0.7,45*91,1);
vals(randsample(1:45*91,3000)) = 0.1;
B = zeros(90); % define B of required size
[ii jj] = ndgrid(1:90); % ii and jj are row and column indices respectively
B(ii>=jj) = vals; % fill in values in column-major order
DTI_mat.one  = B + tril(B,-1).';
vals = unifrnd(0.3,0.7,45*91,1);
vals(randsample(1:45*91,3000)) = 0.1;
B = zeros(90); % define B of required size
[ii jj] = ndgrid(1:90); % ii and jj are row and column indices respectively
B(ii>=jj) = vals; % fill in values in column-major order
DTI_mat.two  = B + tril(B,-1).';
% Make diagonals non-zero
for i = 1:90
    DTI_mat.one(i,i) = min(DTI_mat.one(i,i)+0.5,1);
    DTI_mat.two(i,i) = min(DTI_mat.two(i,i)+0.5,1);
end
% Probit regression sparsity parameter
alpha0_sim =[-2.5,-2.5]; % baseline
% Probit regression beta parameter %
alpha1_sim = [5,5]; % beta
% Sampling of non-zero effective connectivities using standard normal CDF
% probabilities using auxiliary varaible z
z_sim = [];
DTI_vec = [];
gamma_sim =[];
names = {'one','two'};
for i = 1:G
DTI_vec.(names{i})=reshape(DTI_mat.(names{i}),R^2*L,1);
z_sim.(names{i}) =normrnd(alpha0_sim(i)+alpha1_sim(i)*DTI_vec.(names{i}),1,R^2,1);
gamma_sim.(names{i})=(z_sim.(names{i})>0)*1;
end

% Simulate group effective connectivity 
PHI_g=[];
for i = 1:G
    PHI_g.(names{i}) = reshape(gamma_sim.(names{i}).*round(unifrnd(-0.2,0.2,R^2,1),2),R,R);  
end
% Covariance of VAR error term
XI_sim = eye(R);
% Simulate subject level effective connectivities by adding subject-specific
% random matrix to group matrix
% Create a string to store subject names
subjects = cell(n,1);
for i = 1:n
    subjects{i}=strcat('s_',num2str(i));
end

ev = unifrnd(-0.4,0.3,1,90);
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