% Plot steady state objects
gSS = reshape(ggSS,I,2);
g_a = (lla(2) / (lla(1) + lla(2))) * gSS(:,1) ...
    + (lla(1) / (lla(1) + lla(2))) * gSS(:,2);

figure
hold on
f = bar(a / (wSS * (1 - ttau) * zAvg),g_a,'histc');
sh=findall(gcf,'marker','*'); delete(sh);
set(f,'FaceColor',[.03,.24,.8],'EdgeColor',[.8,.24,.03]);
xlim([min(a / (wSS * (1 - ttau) * zAvg)) .8*max(a / (wSS * (1 - ttau) * zAvg))])
xlabel('Net worth to average income','interpreter','latex')
ylabel('Mass of households','interpreter','latex')
title('Steady State Distribution of Assets','interpreter','latex','fontsize',14)
grid on
alpha(0.9)
set(gcf,'color','w')
hold off

% Plot steady state MPCs
tau_MPC = 1;
J = 2;
N_MPC = 21;
dt_MPC = tau_MPC / (N_MPC - 1);

C_stacked = zeros(I*J,N_MPC );
ctilde = zeros(I*J,1);
C_stacked(:,N_MPC) = 0;

B = (1/dt_MPC) * speye(I*J) - A;
c_stacked = reshape(cSS,I*J,1);

for n=N_MPC-1:-1:1
    vec = c_stacked + C_stacked(:,n+1)/dt_MPC;
    C_stacked(:,n) = B\vec;
end

C = reshape(C_stacked(:,1),I,J);
MPC = (C(2:I,:) - C(1:I-1,:)) ./ (aa(2:I,:) - aa(1:I-1,:));

figure
hold on
plot(a(2:I) / (wSS * (1 - ttau) * zAvg),MPC(:,1),'linewidth',1.5,'linestyle','-','color',[.03,.24,.8])
plot(a(2:I) / (wSS * (1 - ttau) * zAvg),MPC(:,2),'linewidth',1.5,'linestyle','-','color',[.8,.24,.03])
plot(a(2:I) / (wSS * (1 - ttau) * zAvg),ones(I-1,1)*MPC(.8*I-1,2),'linewidth',1.5,'linestyle','--','color','k')
xlim([min(a(2:I) / (wSS * (1 - ttau) * zAvg)) .8*max(a / (wSS * (1 - ttau) * zAvg))])
xlabel('Net worth to average income','interpreter','latex')
ylabel('MPC','interpreter','latex')
text(3,.05,'$\leftarrow z=0 $','Fontsize',20,'interpreter','latex','Color',[.8,.24,.03]);
text(1,0.25,'$\leftarrow z=1 $','Fontsize',20,'interpreter','latex','Color',[.03,.24,.8]);
set(gcf,'color','w')
title('MPCs in steady state','interpreter','latex','fontsize',14)
alpha(0.7)
grid on
hold off
drawnow