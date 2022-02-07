%% Settings
clear all;
trials = 1e6; % Number of simulation trails

%% NETWORK SETTINGS
Mx = 16;
My = 16;
M = Mx*My; % Number of reflecting elements

kappa = 1; % Squared amplitude reflection coeff.
q = 2; % Antenna directivity coeff.
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
filename = sprintf('omega_%dx%d.mat',Mx,My);
Omega = cell2mat(struct2cell(load(filename)));

filename = sprintf('amr_%dx%d.mat',Mx,My);
data_amr = cell2mat(struct2cell(load(filename)));

alpha= data_amr(1);
mu = data_amr(2);
r = data_amr(3);

cdf_alpha_mu = @(w) gammainc( mu*(w/r).^(alpha), mu, 'lower' );

% Outage Prob
PTdB = 5:0.25:45; gthdB = 5; gth = 2^(6.0)-1;
for isnr = 1:length(PTdB)
    PT = 10^(PTdB(isnr)/10);
    avgsnr= PT/sigma2;
    
    % = DIRECT LINK SIM =
    OP_drt_sim(isnr) = mean(avgsnr*normsnr_drt < gth);
    
    % = DIRECT LINK ANA =
    cdf_d_drt = @(r) 1/(32*R_icd^3*R_rfl^3) * (32*R_icd^3*r.^3) .* (r <= R_rfl-R_icd)...
            + 1/(32*R_icd^3*R_rfl^3) * (r.^6 - 9*(R_icd^2+R_rfl^2)*r.^4 + 16*(R_icd^3+R_rfl^3)*r.^3 ...
            - 9*(R_icd^2-R_rfl^2)^2*r.^2 + (R_icd-R_rfl)^4*(R_icd^2+4*R_icd*R_rfl+R_rfl^2)).* (r > R_rfl-R_icd).*(r <= R_rfl+R_icd)...
            + (r > R_rfl+R_icd);        
    cdf_snr_drt_D = @(x) 1 - cdf_d_drt( (1/delta)*lambertw( delta*sqrt(avgsnr*chi_drt./x) ) );
    OP_drt_ana(isnr) = cdf_snr_drt_D( gth );
    
    % CASCADED LINK    
    chi_ccd = kappa * (waveLen/(4*pi))^4 * antenna_gain_Tx * antenna_gain_Rx;
    
    OP_rfl_sim(isnr) = mean(avgsnr * chi_ccd * Omega.^2 < gth);
    OP_rfl_ana(isnr) = cdf_alpha_mu( sqrt(gth/avgsnr/chi_ccd) );
end
save('OP_drt_ana.mat','OP_drt_ana');
save('OP_drt_sim.mat','OP_drt_sim');

filename = sprintf('OP_rfl_ana_%dx%d.mat',Mx,My);
save(filename,'OP_rfl_ana');
filename = sprintf('OP_rfl_sim_%dx%d.mat',Mx,My);
save(filename,'OP_rfl_sim');
%%

figure;
semilogy(PTdB,OP_drt_ana,'-k'); hold on;
semilogy(PTdB,OP_rfl_ana,'--k'); hold on;
semilogy(PTdB,OP_drt_sim,'ok'); hold on;
semilogy(PTdB,OP_rfl_sim,'ok'); hold on;

ylabel('Outage Probability','Interpreter','LaTex');
xlabel('$P_{\mathrm{S}}$ [dBm]','Interpreter','LaTex');
legend('Ana. Direct','Sim. Direct','Sim. Cascades');
set(gca, 'LooseInset', get(gca, 'TightInset'));
set(gca,'FontSize',20);
axis([-Inf Inf 1e-2 1]);