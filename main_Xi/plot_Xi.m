clear all;
Xi_icd_sim_2 = cell2mat(struct2cell(load('Xi_icd_sim_q_2.mat')));
Xi_icd_ana_2 = cell2mat(struct2cell(load('Xi_icd_ana_q_2.mat')));

Xi_icd_sim_4 = cell2mat(struct2cell(load('Xi_icd_sim_q_4.mat')));
Xi_icd_ana_4 = cell2mat(struct2cell(load('Xi_icd_ana_q_4.mat')));

Xi_icd_sim_8 = cell2mat(struct2cell(load('Xi_icd_sim_q_8.mat')));
Xi_icd_ana_8 = cell2mat(struct2cell(load('Xi_icd_ana_q_8.mat')));

Xi_icd_sim_16 = cell2mat(struct2cell(load('Xi_icd_sim_q_16.mat')));
Xi_icd_ana_16 = cell2mat(struct2cell(load('Xi_icd_ana_q_16.mat')));

xx = linspace(0,20,25); 
figure;
[y,x] = ecdf(Xi_icd_sim_2);
p1 = plot(x,y,'-r','LineWidth',1.5); hold on;
     plot(xx,Xi_icd_ana_2,'ok');

[y,x] = ecdf(Xi_icd_sim_4);
p2 = plot(x,y,'-b','LineWidth',1.5); hold on;
     plot(xx,Xi_icd_ana_4,'ok');
    
[y,x] = ecdf(Xi_icd_sim_8);
p3 = plot(x,y,'-k','LineWidth',1.5); hold on;
     plot(xx,Xi_icd_ana_8,'ok');

[y,x] = ecdf(Xi_icd_sim_16);
p4 = plot(x,y,'-m','LineWidth',1.5); hold on;
p5 = plot(xx,Xi_icd_ana_16,'ok');

axis([0 20 0.5 1]);

leg = legend([p1 p2 p3 p4 p5],{'q = 2 (Sim)','q = 4 (Sim)','q = 8 (Sim)','q = 16 (Sim)','Ana'});
xlabel('x');
ylabel('CDF');
leg.ItemTokenSize = [10,18];
grid on;
set(gca,'Fontsize',10);
set(gcf,'Position',[100 100 400 200]);
set(gca,'LooseInset',get(gca,'TightInset'));