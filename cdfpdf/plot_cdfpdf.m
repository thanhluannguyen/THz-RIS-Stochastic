% rng('default');
clear all;

Mx = 64;
My = 64;
M = Mx*My;

filename = sprintf('amr_%dx%d.mat',Mx,My);
data_amr = cell2mat(struct2cell(load(filename)));

alpha= data_amr(1);
mu = data_amr(2);
r = data_amr(3);

filename = sprintf('omega_%dx%d.mat',Mx,My);
data_Omega = cell2mat(struct2cell(load(filename)));

%% The CDF
cdf_alpha_mu = @(w) gammainc( mu*(w/r).^(alpha), mu, 'lower' );
% ----------------- CDF
figure(1);
set(gcf,'position',[74.6,380,400,162.4]);
subplot(1,2,1);
[y,xx] = ecdf(data_Omega);
plot(xx,y,'Linewidth',1); hold on;

xxx = logspace(log10(min(xx)),log10(max(xx)),100);
stem(xxx,cdf_alpha_mu(xxx),'markersize',6);
xlabel('$x\times 10^4$','Interpreter','LaTex');
ylabel('CDF','Interpreter','LaTex');
legend('Exact','Approx.');
set(gca,'ytick',[1e-5 1e-4 1e-3 1e-2 1e-1 1]);
set(gca,'xtick', [1.2 1.3 1.4 1.5]*1e4);
set(gca, 'XTickLabel', [1.2 1.3 1.4 1.5]);
set(gca,'XScale','log');
set(gca,'YScale','log');
grid on;
axis([-Inf 1.5e4 1e-5 1]);
%% The PDF
figure(1);
subplot(1,2,2);
histogram(data_Omega','normalization','pdf',...
    'NumBins',150); hold on;

% syms w
% pdf_alpha_mu = @(x) limit( alpha*w^w/(r^(alpha*w)*gamma(w))*x.^(alpha*w-1).* exp(-w*(x/r).^alpha),w,mu );
xxx = xx(1:100:length(xx));
% yyy = double(pdf_alpha_mu(xxx));
plot(xxx,yyy,'linewidth',1.5);
legend('Exact','Approx.');
xlabel('$x\times 10^4$','Interpreter','LaTex');
ylabel('PDF','Interpreter','LaTex');
% set(gca, 'LooseInset', get(gca, 'TightInset'));
axis([-Inf 1.6e4 0 1.2e-3]);
grid on;
set(gca,'xtick', [1.2 1.3 1.4 1.5 1.6]*1e4);
set(gca, 'XTickLabel', [1.2 1.3 1.4 1.5 1.6]);
set(gca,'ytick',[0 0.25 0.5 0.75 1]*10^(-3));