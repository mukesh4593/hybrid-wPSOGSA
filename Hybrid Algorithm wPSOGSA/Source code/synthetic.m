% close all;clear all;clc
% code for generating synthetic MT data
% Written/modified by Mukesh Mukesh, Kuldeep Sarkar and Upendra K. Singh in 2019
%% Input
% resistivities: true resistivity of various layer
% thicknesses: thickness of various layers
% frequencies: frequency range for generating MT synthetic data
% nnn: No. of observations or the number of frequencies in frequency range

%% Output
% Observe apparent resistivity
% Observe apparent phase

frequencies=[5.88,4.55,3.03,2.22,1.47,1.14,0.741,0.552,0.377,0.141,0.093,0.0709,0.0466,0.0348,0.0234,0.0175,0.0117,0.00875,0.00587,0.00439,0.00284,0.00216,0.00145,0.00106,0.000674,0.000486];
resistivities = [30000 5000 1000];% True model taken from shaw & shalivahan (2007)
thicknesses = [15000 18000];
nnn=length(frequencies);
for k = 1:nnn
    frequency = frequencies(k);
    mu = 4*pi*1E-7; %Magnetic Permeability (H/m)  
w = 2 * pi .* frequency; %Angular Frequency (Radians);
n=length(resistivities); %Number of Layers
impedances = zeros(n,1);

% Calculate basement impedance  
Zn = sqrt(sqrt(-1)*w*mu*resistivities(n)); 
impedances(n) = Zn; 
%Iterate through layers starting from layer j=n-1 (i.e. the layer above the basement)        
for j = n-1:-1:1
    resistivity = resistivities(j);
    thickness = thicknesses(j);
    dj = sqrt(sqrt(-1)* (w * mu * (1/resistivity)));
    wj = dj * resistivity;
       ej = exp(-2*thickness*dj);                     
      belowImpedance = impedances(j + 1);
    rj = (wj - belowImpedance)/(wj + belowImpedance); 
    re = rj*ej; 
    Zj = wj * ((1 - re)/(1 + re));
    impedances(j) = Zj;               
end
% Compute apparent resistivity from top layer impedance
Z= impedances(1);
absZ = abs(Z);
apparentResistivity = (absZ * absZ)/(mu * w);

r_cal(1,k) =frequency;
r_cal(2,k) = apparentResistivity;
re=real(Z);
img=imag(Z);
x=atan(img./re);
r_cal(3,k)=x;%.*(180/pi);
end
obs_data11=r_cal';
obs_data11=abs(obs_data11);%%%% Absolute value 
obs_data11(:,1)=frequencies;

figure
loglog(frequencies,obs_data11(:,2),'*')
r_cal=r_cal';
xlabel('Frequency (Hz)')
ylabel('Apparent Resistivity (ohm-m)') 
title('Synthetic Model')
%%%%%%Phase Calculation
phase=obs_data11(:,3).*(180/pi);%%Phase in degree
figure
semilogx(frequencies,phase)
xlabel('Frequency (Hz)')
ylabel('Phase') 
%%%%%%%save file in (.dat) format
save obs_data11.dat obs_data11 -ascii
