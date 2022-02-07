figure;
plot(xx,yy,'Linewidth',1); hold on;
plot(xxx,yyy,'--','Linewidth',1.5);
xlabel('$x$','Interpreter','LaTex');
ylabel('CDF','Interpreter','LaTex');
legend('Exact','Approximation','Orientation','Horizontal');
set(gca,'XScale','log');
set(gca,'YScale','log');
grid on;

switch M
    case 16*16
        axis([-Inf 1000 1e-5 1]);
    case 32*32
        axis([-Inf 4000 1e-5 1]);
    case 64*64
        axis([-Inf 1.5e4 1e-5 1]);
    case 96*96
        axis([-Inf 1.7e4 1e-5 1]);
end