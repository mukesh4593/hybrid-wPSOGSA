% Hybrid wPSOGSA algorithm
clear all
close all
clc
data = load('obs_data11.dat');%Observed Data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dataFrequencies =data(:,1);% 1./period;
r_obs= data(:,2);% observed apparent resistivity
p_obs= data(:,3);%observed apparent phase
N = 50;     % Size of the swarm " no agents/particles "
Max_Iteration =1000;    % Maximum number of "iterations" % MaxIt=Max_Iteration;  % PSO
dim=5;% No. of layer parameters.
run=10;% Number of computations/Models
%% Search range
down=[5000  1000  50    5000   10000];
up=[50000	10000 5000	25000  25000];
%%

%%%% wPSOGSA
for k=1:run
  k
[gBestScore,gBest,GlobalBestCost,gbest1,r_calPG,p_calPG]= PSOGSA(N, Max_Iteration,dataFrequencies,r_obs,down,up,dim,p_obs);
gbest_run(k,:)=gBest;
gBestScore_run(k,:)=gBestScore;
GlobalBestCost_run(k,:)=GlobalBestCost;
gbest1_run(k).m=gbest1;%%% PArameter for All models
r_cal_PG(:,k)=r_calPG;
p_cal_PG(:,k)=p_calPG;
end
[gbscore,indexPG]=min(gBestScore_run); % finding index for least misfit error
gbestmodel=gbest_run(indexPG,:);% best model based on misfit error
globalbestcost=GlobalBestCost_run(indexPG,:); % least misfit error
gbest11=gbest1_run(indexPG).m; % model saved at each iteration
r_cal11=r_cal_PG(:,indexPG); % calculated apparent resistivity based on misfit error
p_cal11=p_cal_PG(:,indexPG); % calculated apparent phase based on misfit error
%%%END
%% Plot
figure
loglog(dataFrequencies,r_obs,'r*',dataFrequencies,r_cal11,'g')
xlabel('Data Frequency (Hz)');ylabel('Apparent Resistivity (ohm-m)');
legend('Observed Data','Calculated Data');
title('Apparent Resistivity Curve of wPSOGSA')

% %%%%Plot of Convergence Rate
figure
semilogy(globalbestcost,'-r');%wPSOGSA
xlabel('No. of Iteration');ylabel('Misfit');set(gca, 'YScale', 'log')
legend('\fontsize{12} wPSOGSA');
