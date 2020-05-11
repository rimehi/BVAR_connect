% This script is used to generate simulated data with different scenarios 
% and assess performance of the VB model using a set of classification
% metrics. For the sake of demonstration, we assume that there are two
% groups of population.

%% First simulate data using one of the four available options.
seed = 2; % Define a random number for the seed.
data = sim_simulation(seed); % Code for R = 10, L = 1 
%data = sim_simulation_lag_2(seed); % Code for R = 10 , L = 2
% data = sim_simulation_large_30(seed); %Code for R = 30, L =1
% data = sim_simulation_large_90(seed); %Code for R = 90, L =1
% Input correct number of lags and ROIs here.
L = 1; % Number of Lags
G = 2; % Number of Groups
n = length(data.eta);
R = 10; % Number of Regions of Interest
% R = 30;
% R = 90;
T = size(data.X,1);
% Prior to using the simulated data, make sure that data set is stationary
for i = 1:length(data.eta)
    min_val = min(min(data.X(:,:,i)));
    max_val = max(max(data.X(:,:,i)));
    diff(i) = max_val-min_val;
end
if max(diff) > 20
    sprintf('Data is non-stationary, consider using a different seed')
else
    sprintf('Data is stationary, proceed to next step with simulated data')
end
% If the data is not considered to be stationary based on this simpl
% heuristics, try a different number for the seed before proceeding to
% the next step.
%% Create .mat files to be used as inputs the BVAR connect software.
X = data.X;
ROI_names = cell(1,R);
for i = 1:R
    ROI_names{i} = strcat('ROI_',num2str(i));
end
eta = data.eta;
save('fmri_dat.mat','X','ROI_names','L','G','eta')
DTI_vec =cell(1,2); % Matrix defining structural connectivity
DTI_vec{1} = data.DTI_vec.one;
DTI_vec{2} = data.DTI_vec.two;
save('DTI_vec.mat','DTI_vec')
%% Run the BVAR GUI
% Set q to large value
% Use default ICAR prior setting.
% Once code has successfully run, load the output.mat file
load output.mat 
% Evaluate model performance
eval_stat(data,out)

