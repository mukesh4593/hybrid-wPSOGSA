% Code of function for calculating posterior PDF
% Developed by Mukesh Mukesh, kuldeep Sarkar and Upendra K. Singh in 2019

% dim: number of layer parameters
% model:total number of runs/ total generated models after inversion
% percentage:desired confidence interval in percentage
% r_obs: observed apparent resistivity data
% GlobalBestCost_run:misfit error at each iteration for all runs
% gbest_run:output layer parametsrs for all runs
% r_cal_PG:calculated apparent resistivity for all runs
% threshold: threshold misfit error 



function[ac68_pos,bc68_pos]=post(dim,model,percentage,r_obs,GlobalBestCost_run,gbest_run,r_cal_PG,threshold);
m=length(r_obs);
x = GlobalBestCost_run(:,1000); % Misfit error of last iteration for all models
z = find(abs(x)<threshold); % Input the threshold error 
lbest=gbest_run; % output layer parameters for all models
r_cal_c=r_cal11; % Calculated apparent resistivity for all models

  k=0;
for k=1:dim
a(k)=mean(lbest(:,k));
b(:,k)=(lbest(:,k)-a(k));
end
b=b.^2;
  k=0;
for k=1:model
e(:,k)= r_obs-r_cal_c(:,k); 
ee(k)=mean(e(:,k));
eee(:,k)=(e(:,k)-ee(k));
end
eee=eee.^2;

k=0;
for k=1:dim
pd1(:,k)=fitdist(lbest(:,k),'Normal');
y1(:,k)=pdf(pd1(:,k),lbest(:,k));
end

k=0;
for k=1:model
phi=mean(e(:,k).^2);
nor=norm(e(:,k),1);
yyy(k)=((1/(sqrt(2*pi)*phi))^m/2).*exp((-(nor)^2)/(2*m*phi));
end
%%%%%%%%% Likelihood
k=0;
for k=1:model
    lik(k,:)=yyy(k).*y1(k,:);
    slik(k)=sum(lik(k,:));
end
k=0;j=0;
for k=1:model
    for j=1:dim
        pos(k,j)=lik(k,j)./slik(k);
    end
end
%%%%%%% Threshold Data for 68.27% confidence interval
j=0;
   for j=1:length(z)
    kt(j,:)=lbest(z(j),:);
   pos1(j,:)=pos(z(j),:);
   end
   %% Calculation of minimum and maximum confidence interval
meanh=mean(kt);
sdh=std(kt);
ww=(sdh./100);
sd68=ww.*percentage;
mean68p=meanh+sd68;
mean68m=meanh-sd68;
   %%
%    Plot
   j=0;
   figure
for j=1:dim
    subplot(2,dim,j)
    plot(kt(:,j),pos1(:,j),'.')
end
%% Posterior PDF calculation for 68.27% confidence interval
k=0;j=0;
for k=1:dim
    for j=1:length(kt)
        if kt(j,k)>=mean68m(1,k) & kt(j,k)<=mean68p(1,k)
            y85(j,k)=kt(j,k);
        else 
            y85(j,k)=0;
        end
    end
end
k=0;j=0;
for k=1:dim
    for j=1:length(kt)
        if y85(j,k)~=0
            y850(j,k)=pos1(j,k);
        else 
            y850(j,k)=0;
        end
    end
end


k=0;
j=0;
for k=1:dim
    xxx=[];
    xx=[];
    yy=[];
    xxx=find(y85(:,k));
    for j=1:length(xxx)
       z1= xxx(j,1);
    xx(j,1)=y85(z1,k);
    yy(j,1)=y850(z1,k);
    end
    mmmm(k).aa(:,1)=xx;
    mmmm(k).aa(:,2)=yy;
    subplot(2,dim,k+dim)
    plot(xx,yy,'.')
    ac68_pos(k)=mean(xx); % mean model
    bc68_pos(k)=std(xx);  % uncertainity
end
