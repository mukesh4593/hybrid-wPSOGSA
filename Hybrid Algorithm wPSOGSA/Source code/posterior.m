%%% Run/Execute after data inversion

% Threshold is the threshold error while estimating posterior PDF
% r_cal_PG is the calculated apparent resistivity for given number of run/models
% gbest_run is the output layer parameters for required runs
% GlobalBestCost_run is the misfit for required runs
% r_obs is the observed apparent resistivity named as data
dim=5; %%%%%% No. of parameter
run=10; %%% No. of model or run
percentage=68.27; %%%confident interval
threshold=0.0001;% threshold error
%%load files: r_obs, GlobalBestCost_run, gbest_run, r_cal_PG

[ac68_pos,bc68_pos]=post(dim,run,percentage,r_obs,GlobalBestCost_run,gbest_run,r_cal_PG,threshold);
