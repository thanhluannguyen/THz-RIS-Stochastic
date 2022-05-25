%% Settings
clear all;
trials = 1e6; % Number of simulation trails

%% NETWORK SETTINGS
Mx = 16;
My = 16;
M = Mx*My; % Number of reflecting elements

kappa = 1; % Squared amplitude reflection coeff.
q = 16; % Antenna directivity coeff.
epsilon = 1; % Antenna efficiency

% Carrier frequency (in GHz)
fc = 300*10^9; % Hz = 300 GHz
absLossCoeff = 0.0033; % m^{-1} Absorbtion Loss Coefficient Measured at 300 GHz
delta = absLossCoeff/2;
waveLen = physconst('LightSpeed')/fc;

% Network area
Radius = 1; % Unit
x_area_min = -Radius;
x_area_max = Radius;
y_area_min = -Radius;
y_area_max = Radius;
z_area_min = -Radius;
z_area_max = Radius;

% Location of the m-th reflecing element
x_ris = 0;
y_ris = 0;
z_ris = 0;

% Maximum distances to the origin
R_icd = y_area_max/2;
R_rfl = y_area_max;

vphi_icd = 2*pi*rand(1,trials);
cospsi_icd = -1 + 2*rand(1,trials);
psi_icd = acos(cospsi_icd);
d_icd = R_icd * rand(1,trials).^(1/3);

x_Tx = d_icd.*cos(vphi_icd).*sin(psi_icd);
y_Tx = d_icd.*sin(vphi_icd).*sin(psi_icd);
z_Tx = d_icd.*cos(psi_icd);

vphi_rfl = 2*pi*rand(1,trials);
cospsi_rfl = -1 + 2*rand(1,trials);
psi_rfl = acos(cospsi_rfl);
d_rfl = R_rfl * rand(1,trials).^(1/3);

x_Rx = d_rfl.*cos(vphi_rfl).*sin(psi_rfl);
y_Rx = d_rfl.*sin(vphi_rfl).*sin(psi_rfl);
z_Rx = d_rfl.*cos(psi_rfl);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
antenna_gain_Tx = db2pow(0); % Antenna gain of Tx
antenna_gain_Rx = db2pow(0); % Antenna gain at Rx
antenna_gain_icd = @(psi,phi) 2*epsilon*(q+1)*(cos(psi).^q).*(psi < pi/2).*(psi > 0); % Antenna gain of each reflecting element
antenna_gain_rfl = @(psi,phi) 2*epsilon*(q+1)*(cos(psi).^q).*(psi < pi/2).*(psi > 0); % Antenna gain of each reflecting element

% Bandwidth
BW = 100*1e6; % 100 MHz 

% Noise figure (in dB)
noiseFiguredB = 6.5; % dB THz, from 280 GHz - 330 GHz

% Compute the noise power in dBm
sigma2dBm = -174 + 10*log10(BW) + noiseFiguredB; % dBm
sigma2 = db2pow(sigma2dBm);

%% Normalized SNR over the Direct Tx-Rx link
D_drt = sqrt( (x_Rx-x_Tx).^2+(y_Rx-y_Tx).^2+(z_Rx-z_Tx).^2 );

PL_drt= antenna_gain_Tx .* antenna_gain_Rx * (waveLen/(4*pi))^2 ...
        .* exp( -delta*D_drt )./D_drt.^2;
    
h_drt = sqrt(PL_drt) .* exp( -1i*2*pi/waveLen*D_drt ); % Assuming Perfect Phase-Shift
chi_drt = (waveLen/(4*pi))^2 * antenna_gain_Tx * antenna_gain_Rx;

normsnr_drt = abs(h_drt).^2;

% Cascaded Links
tmp = sqrt(2*epsilon*(q+1));

a_xi = @(xi) min(tmp,R_icd*exp(delta*R_icd)*xi);
b_xi = @(xi) lambertw(delta*a_xi(xi)./xi);

Xi_icd_sim = sqrt(antenna_gain_icd(psi_icd,vphi_icd)).*exp(-delta*d_icd)./d_icd;
    
Xi_icd_ana = @(xi) 1/2 + 3/(2*delta^3*R_icd^3)*(xi/tmp/delta).^(2/q)...
    * (-1)^(-2/q+1)/(2/q)^(2/q+3) * gammainc(-2/q*b_xi(xi),2/q+3,'lower')*gamma(2/q+3)...
    + ( 1 - (b_xi(xi)/(delta*R_icd)).^3 )/2.*(xi > exp(-delta*R_icd)/R_icd*tmp);
    
xx = linspace(0,20,25); 
fXi_icd_ana = zeros(size(xx));
for ix = 1:length(xx)
    fXi_icd_ana(ix) = Xi_icd_ana(xx(ix));
end

filename = sprintf('Xi_icd_sim_q_%d.mat',q);
save(filename,'Xi_icd_sim');
filename = sprintf('Xi_icd_ana_q_%d.mat',q);
save(filename,'fXi_icd_ana');

figure;
[y,x] = ecdf(Xi_icd_sim);
plot(x,y,'-b','LineWidth',1.5); hold on;
plot(xx,fXi_icd_ana,'ok');
axis([0 20 0.5 1]);
xlabel('x');
ylabel('CDF');