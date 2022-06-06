% Function of calculating cost function 
% Developed/modified by Mukesh Mukesh, kuldeep Sarkar and Upendra K. Singh for geophysical data for academic purpose in 2019.
% % Input
% ginv:inverted data
% gobs: observed data

% % Output
% Fitness
function fitness=RMS_1(ginv,gobs)
sum1=0;
sum2=0;
for k=1:length(gobs)
    sum1=sum1+abs(ginv(k)-gobs(k));
    sum2=sum2+abs(ginv(k)+gobs(k));
end
fitness=2*sum1/(sum1+sum2);
% 
% e=gobs-ginv;
% 
% fitness=norm(e);
