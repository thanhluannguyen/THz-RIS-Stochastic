%% Settings
clear all;
trials = 5e6; % Number of simulation trails

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

pdf_d_drt = @(r) 3/(16*R_icd^3*R_rfl^3) * (16*R_icd^3*r.^2) .* (r <= R_rfl-R_icd)...
        + 3/(16*R_icd^3*R_rfl^3) * (r.^5 - 6*(R_icd^2+R_rfl^2)*r.^3 + 8*(R_icd^3+R_rfl^3)*r.^2 ...
            - 3*(R_icd^2-R_rfl^2)^2*r).* (r > R_rfl-R_icd).*(r <= R_rfl+R_icd);
cdf_d_drt = @(r) 1/(32*R_icd^3*R_rfl^3) * (32*R_icd^3*r.^3) .* (r <= R_rfl-R_icd)...
        + 1/(32*R_icd^3*R_rfl^3) * (r.^6 - 9*(R_icd^2+R_rfl^2)*r.^4 + 16*(R_icd^3+R_rfl^3)*r.^3 ...
        - 9*(R_icd^2-R_rfl^2)^2*r.^2 + (R_icd-R_rfl)^4*(R_icd^2+4*R_icd*R_rfl+R_rfl^2)).* (r > R_rfl-R_icd).*(r <= R_rfl+R_icd)...
        + (r > R_rfl+R_icd);

cdf_normsnr_drt = @(x) 1 - cdf_d_drt( (1/delta)*lambertw( delta*sqrt(chi_drt./x) ) );

pdf_normsnr_drt = @(x) lambertw( delta*sqrt(chi_drt./x) )...
    ./ ( 1+lambertw( delta*sqrt(chi_drt./x) ) )./(2*x*delta)...
        .* pdf_d_drt( (1/delta)*lambertw( delta*sqrt(chi_drt./x) ) );

%%
figure;
[y,x] = ecdf(normsnr_drt);
plo_1 = plot(x,y,'linewidth',2); hold on;
plot(x,cdf_normsnr_drt(x),'--y','linewidth',1.5); hold on;

[y,x] = ecdf(db2pow(5)*normsnr_drt);
plo_2 = plot(x,y,'linewidth',2); hold on;
plot(x,cdf_normsnr_drt(x/db2pow(5)),'--y','linewidth',1.5); hold on;

[y,x] = ecdf(db2pow(10)*normsnr_drt);
plo_3 = plot(x,y,'linewidth',2); hold on;
plot(x,cdf_normsnr_drt(x/db2pow(10)),'--y','linewidth',1.5); hold on;

[y,x] = ecdf(db2pow(15)*normsnr_drt);
plo_4 = plot(x,y,'linewidth',2); hold on;
plot(x,cdf_normsnr_drt(x/db2pow(15)),'--y','linewidth',1.5); hold on;

[y,x] = ecdf(db2pow(20)*normsnr_drt);
plo_5 = plot(x,y,'linewidth',2); hold on;
plo_6 = plot(x,cdf_normsnr_drt(x/db2pow(20)),'--y','linewidth',1.5); hold on;

set(gca,'XScale','log');
set(gca,'YScale','log');
leg = legend([plo_1 plo_2 plo_3 plo_4 plo_5 plo_6],...
    'SNR = 0 dB (Sim)',...
    'SNR = 5 dB (Sim)',...
    'SNR = 10 dB (Sim)',...
    'SNR = 15 dB (Sim)',...
    'SNR = 20 dB (Sim)',...
    'Analytical','NumColumns',2);
leg.ItemTokenSize = [10,18];
grid on;
set(gca,'Fontsize',10);
set(gcf,'Position',[100 100 400 200]);
axis([10^(-9) 10^(-6) 10^(-5.5) 1]);
set(gca,'LooseInset',get(gca,'TightInset'));
xlabel('$x$','Interpreter','LaTex');
ylabel('CDF','Interpreter','LaTex');

%%
figure;
% histogram(normsnr_drt,'normalization','pdf','NumBins',10); hold on;

scale = 1e3;
% [y,x] = ksdensity(normsnr_drt,'Numpoints',5e6); hold on;

set(gca,'XScale','log');

xx = logspace(-9,-6,1e4);
yy = pdf_normsnr_drt(xx);
histogram(normsnr_drt,'normalization','pdf','NumBins',3e5); hold on;
% plot(x,y); hold on;
plot(xx,yy,'-r','linewidth',1.5); hold on;
set(gca,'XScale','log');
leg = legend('Simulation','Analytical');
axis([10^(-8.8) 10^(-7) -Inf 11*10^7]);

leg.ItemTokenSize = [10,18];
grid on;
set(gca,'Fontsize',10);
set(gcf,'Position',[100 100 400 200]);
set(gca,'LooseInset',get(gca,'TightInset'));

xlabel('$x$','Interpreter','LaTex');
ylabel('PDF','Interpreter','LaTex');