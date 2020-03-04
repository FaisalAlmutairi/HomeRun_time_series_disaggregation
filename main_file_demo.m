%% main file of the code (demo)
%  This code performs an example of the experiments in the paper titled:
%  "HomeRun: Scalable Sparse-Spectrum Reconstruction of Aggregated Historical Data"
%
%  Here we compare the proposed method HomeRun (HR) to the baselines: 
%  H-Fuse [Liu et al. SDM'17], and Least Squares under various scenarios 
%  using measles counts data (https://www.tycho.pitt.edu).
%  
%  This code produces Figure 9 in the paper.
%
%  Faisal Almutairi (almut012@umn.edu) , and Fan Yang (fay28@pitt.edu), Jan. 2018.


clear all; close all; clc; load NYC_measles_counts.mat;

x = NYC_measles_counts';
N = length(x);

% create the smoothness matrix, named H
h = [1 -1];
c = [h(1); zeros(N-2,1)];
r = zeros(1,N);
r(1:2) = h.';
H = toeplitz(c,r);

% Report Duration (RD) and Shift limits
MinRD = 2;
MaxRD = 52;
MinShift = 1;
MaxShift = 26;
%% Test HomeRun and baselines on various aggregation scenarios        
for RD = MinRD:10:MaxRD
    for Shift = MinShift:2:MaxShift
        fprintf('RD = %d, Shift = %d \n', RD, Shift);
        
        % creat an aggregation (observation) matrix, named O
        Overlap = RD-Shift;
        fit = 1;
        O = create_obs_matrix(N,RD,Overlap,fit);

        
        % create aggregated reports from the time series x, named y
        y = O*x';

        % disaggregate using the baselines H-Fuse, and Least Squares (LS) 
        x_LS  = (pinv(O)*y).';
        x_Hfuse = H_Fuse(O, y, H);


        % disaggregate using the proposed method HomeRun (HR)
        [s_HR, cost] = (HomeRun(O, H, y, 1, 1));
        % compute the inverse Discrete Cosine Transform (function is written by Andriy Myronenko)
        x_HR = mirt_idctn(full(s_HR)); 
        
        % computing the errors 
        rmse_LS(RD-1,Shift) = sqrt(mean((x-x_LS).^2));
        rmse_Hfuse(RD-1,Shift) = sqrt(mean((x-x_Hfuse).^2));
        rmse_HR(RD-1,Shift) = sqrt(mean((x-x_HR').^2));
        
        % computing the ration of error difference: 
        % + means HomeRun is better, - means baseline is better
        Compare_HR_LS(RD-1,Shift) = (rmse_LS(RD-1,Shift)-rmse_HR(RD-1,Shift))...
                                    /max(rmse_LS(RD-1,Shift),rmse_HR(RD-1,Shift));
        Compare_HR_Hfuse(RD-1,Shift) = (rmse_Hfuse(RD-1,Shift)-rmse_HR(RD-1,Shift))...
                                    /max(rmse_Hfuse(RD-1,Shift),rmse_HR(RD-1,Shift));
    end  
end  
%% generate plots (Fig. 9 in the paper)
figure
ConstrName = 'HomeRun vs. H-Fuse';
plot_error_ratio(Compare_HR_Hfuse,ConstrName);
grid on
figure
ConstrName = 'HomeRun vs. LS';
plot_error_ratio(Compare_HR_LS,ConstrName);
grid on