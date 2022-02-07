% rng('default');
clear all;
close all;

Mx = 64;
My = 64;
M = Mx*My;

filename = sprintf('omega_%dx%d.mat',Mx,My);
Data = cell2mat(struct2cell(load(filename)));

switch M
    case 16*16
        E1 = 7587876969173823/8796093022208;
        E2 = 1631105307663041/2147483648;
        E4 = 5215395313914319/8192;
	case 32*32
        E1 = 948460287810735/274877906944;
        E2 = 6425046437849897/536870912;
        E4 = 2342660441287465/16;
    case 64*64
        E1 = 7587237814760707/549755813888;
        E2 = 6399440580156255/33554432;
        E4 = 36566547396404704;
    case 96*96
        E1 = 8535660068500679/549755813888;
        E2 = 4049532499688133/16777216;
        E4 = 58564876858984312;
end        
C = [E1^2/(E2-E1^2), E2^2/(E4-E2^2), E1];

% Choosing the start points for the solver "fsolve"
% nrmv: normalization vector, reduce computational complexity
switch M
    case 16*16
        nrmv= [ 4970597207513969/140737488355328, 8366278490360667/9007199254740992, 8072794775984711/140737488355328];
        x0 = [ 6990272365865527/70368744177664, 8416214275548159/4611686018427387904, 3783882971820837/9223372036854775808];
    case 32*32
        x0 = [ 16093995955225/70368744177664, 859960919928683/2251799813685248, 1993235507137207/4503599627370496];
        nrmv= [ 5302894522532935/8796093022208, 1, 3019669095777541/17592186044416];
    case 64*64
        x0 = [ 7664658955855355/140737488355328, 1076793298464919/1125899906842624, 2782720393985599/9007199254740992];
        nrmv=[ 8272302692080273/8796093022208, 1, 187709217667025/68719476736];
    case 96*96
        x0 = [ 50, 2727911556966581/2305843009213693952, 7360583020339433/2305843009213693952];
        nrmv= [ 1586523710967105/137438953472, 1, 4537732047708801/2199023255552];
    otherwise
        x0 = [randi([10 100]) rand/randi([10 100]) randi([10 100])];
        nrmv = [rand*1e5 rand*100 rand*1e4];
end
%% DO NOT MODIFY THIS PART
fprintf('Start: %.1f, %.4f, %.4f \n',x0(1),x0(2),x0(3))
fun = @(x) alphamuestm(abs(nrmv.*x),C); % Estimators of alpha-mu distribution
options = optimoptions(@fsolve,'Algorithm','trust-region-dogleg',...
    'MaxFunctionEvaluations',1e6,'MaxIterations',1e4,...
    'TolFun',1e-18,'TolX',1e-18);
[x,fval] = fsolve(fun,x0,options);

mu = nrmv(1)*abs(x(1));
alpha = 1/( nrmv(2)*abs(x(2)) );
r = nrmv(3)*abs(x(3));
%%
amr = [alpha,mu,r];
filename = sprintf('amr_%dx%d.mat',Mx,My);
save(filename,'amr');
%% PLOT THE RESULTS (CDF)
cdf_alpha_mu = @(w) gammainc( mu*(w/r).^(alpha), mu, 'lower' );
% ----------------- CDF
[yy,xx] = ecdf(Data);

nzy = yy(yy > 0); % non-zero-y
nzx = xx(yy > 0); % non-zero-x

xxx = nzx(1:50:length(nzx));
yyy = cdf_alpha_mu(xxx);
yyyy= nzy(1:50:length(nzx));

RMSE = mean( abs(yyyy-yyy)./yyyy )
% if (RMSE > 1e-3) || (isnan(RMSE))
%     main_fit2alphamu;
% else
    %%
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
% end