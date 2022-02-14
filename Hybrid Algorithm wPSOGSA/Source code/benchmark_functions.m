function [fit,r_cal,p_cal]=benchmark_functions(mm,dataFrequencies,r_obs,p_obs)

nnn =length(r_obs);
for k = 1:nnn
a=round(length(mm)/2);
[r_cal(k,1),p_cal(k,1)] = forward(mm(1:a), mm(a+1:length(mm)),dataFrequencies(k));

end

fit=RMS_1(r_cal,r_obs); % Cost function
end
