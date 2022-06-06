% Code for calculating posteior PDF
% Written by Mukesh Mukesh, kuldeep Sarkar and Upendra K. Singh in 2019
% Run/Execute after data inversion

% dim: number of layer parameters
% model:total number of runs/ total generated models after inversion
% percentage:desired confidence interval in percentage
% r_obs: observed apparent resistivity data
% GlobalBestCost_run:misfit error at each iteration for all runs
% gbest_run:output layer parametsrs for all runs
% r_cal_PG:calculated apparent resistivity for all run/models
% threshold: threshold misfit error while estimating posterior PDF
%ac68_pos: mean model for posterior PDF
%bc68_pos: uncertainty for posterior PDF


dim=5; %%%%%% No. of parameter
run=10; %%% No. of model or run
percentage=68.27; %%%confident interval
threshold=0.0001;% threshold error
%%load files: r_obs, GlobalBestCost_run, gbest_run, r_cal_PG

[ac68_pos,bc68_pos]=post(dim,run,percentage,r_obs,GlobalBestCost_run,gbest_run,r_cal_PG,threshold);
