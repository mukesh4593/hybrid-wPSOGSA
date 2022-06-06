% Code for forward modelling of 1D MT inversion
% Developer-This code was Generated by SeyedAli Mirjalili, 2011, A New Hybrid PSOGSA for mathematical benchmark function 
% modified by Mukesh Mukesh, kuldeep Sarkar and Upendra K. Singh for geophysical data for academic purpose in 2019.

% resistivities: true resistivity of various layer
% thicknesses: thickness of various layers
% frequency: frequency range for generating MT synthetic data

function  [apparentResistivity,phase] = forward(resistivities, thicknesses,frequency)
mu = 4*pi*1E-7; %Magnetic Permeability (H/m)  
w = 2 * pi * frequency; %Angular Frequency (Radians);
n=length(resistivities); %Number of Layers
impedances = zeros(n,1); 
Zn = sqrt(sqrt(-1)*w*mu*resistivities(n)); % Impedance calculation
impedances(n) = Zn; 
for j = n-1:-1:1
    dj = sqrt(sqrt(-1)* (w * mu * (1/resistivities(j))));
    wj = dj * resistivities(j);
    ej = exp(-2*thicknesses(j)*dj);                     
    belowImpedance = impedances(j + 1);
    rj = (wj - belowImpedance)/(wj + belowImpedance); 
    re = rj*ej; 
    Zj = wj * ((1 - re)/(1 + re));
    impedances(j) = Zj;               
end
Z = impedances(1);
apparentResistivity = (abs(Z) * abs(Z))/(mu * w); % Apparent resistivity calculation
phase=atan(imag(Z)./real(Z)); % Apparent phase calculation
end    