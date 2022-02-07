blue1=[0, 0.4470, 0.7410];
orange1=[0.8500, 0.3250, 0.0980];
yellow1=[0.9290, 0.6940, 0.1250];
purple1=[0.4940, 0.1840, 0.5560];
green1=[0.4660, 0.6740, 0.1880];
cyan1=[0.3010, 0.7450, 0.9330];
red1=[0.6350, 0.0780, 0.1840];

OP_drt_ana = cell2mat(struct2cell(load('OP_drt_ana.mat')));
OP_drt_sim = cell2mat(struct2cell(load('OP_drt_sim.mat')));

OP_rfl_ana_16 = cell2mat(struct2cell(load('OP_rfl_ana_16x16.mat')));
OP_rfl_ana_32 = cell2mat(struct2cell(load('OP_rfl_ana_32x32.mat')));
OP_rfl_ana_64 = cell2mat(struct2cell(load('OP_rfl_ana_64x64.mat')));
OP_rfl_ana_96 = cell2mat(struct2cell(load('OP_rfl_ana_96x96.mat')));

OP_rfl_sim_16 = cell2mat(struct2cell(load('OP_rfl_sim_16x16.mat')));
OP_rfl_sim_32 = cell2mat(struct2cell(load('OP_rfl_sim_32x32.mat')));
OP_rfl_sim_64 = cell2mat(struct2cell(load('OP_rfl_sim_64x64.mat')));
OP_rfl_sim_96 = cell2mat(struct2cell(load('OP_rfl_sim_96x96.mat')));

PTdB = 5:0.25:45; 
pnt = [1:2:(length(PTdB)-4) (length(PTdB)-3):1:length(PTdB)];
figure;
semilogy(0,0,'--k'); hold on;
semilogy(0,0,'-k');
semilogy(0,0,'ok');

set(gcf,'position',[74.6,400,400,160]);
semilogy(PTdB,OP_drt_ana,'--k','Linewidth',1); hold on;
semilogy(PTdB(pnt),OP_drt_sim(pnt),'ok','MarkerSize',5); hold on;

semilogy(PTdB,OP_rfl_ana_16,'-k','Linewidth',1); hold on;

semilogy(PTdB(pnt),OP_rfl_sim_16(pnt),'ok','MarkerSize',5); hold on;
semilogy(PTdB(pnt),OP_rfl_sim_32(pnt),'ok','MarkerSize',5); hold on;
semilogy(PTdB(pnt),OP_rfl_sim_64(pnt),'ok','MarkerSize',5); hold on;
semilogy(PTdB(pnt),OP_rfl_sim_96(pnt),'ok','MarkerSize',5); hold on;

semilogy(PTdB,OP_rfl_ana_32,'-k','Linewidth',1); hold on;
semilogy(PTdB,OP_rfl_ana_64,'-k','Linewidth',1); hold on;
semilogy(PTdB,OP_rfl_ana_96,'-k','Linewidth',1); hold on;

ylabel('Outage Probability','Interpreter','LaTex');
xlabel('$P_{\mathrm{S}}$ [dBm]','Interpreter','LaTex');
legend('Direct link (ana.)','Cascaded link (ana.)','Sim.');
set(gca, 'LooseInset', get(gca, 'TightInset'));
axis([5 40 1e-3 1]);