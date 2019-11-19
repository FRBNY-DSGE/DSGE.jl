figure
subplot(2,2,1)
hold on
plot(vTime,100 * vAggregateTFP_full,'linewidth',1.5,'linestyle','-','color',[.03,.24,.8])
plot(vTime,100 * vAggregateTFP_reduced,'linewidth',1.5,'linestyle','--','color',[.8,.24,.03])
set(gcf,'color','w')
xlim([vTime(1) vTime(end)])
grid on
title('TFP','interpreter','latex','fontsize',14)
ylabel('$\%$ deviation from s.s.','interpreter','latex')
legend('Full model','Reduced model','location','northeast')
hold off

subplot(2,2,2)
hold on
plot(vTime,100 * vAggregateOutput_full,'linewidth',1.5,'linestyle','-','color',[.03,.24,.8])
plot(vTime,100 * vAggregateOutput_reduced,'linewidth',1.5,'linestyle','--','color',[.8,.24,.03])
set(gcf,'color','w')
xlim([vTime(1) vTime(end)])
grid on
title('Output','interpreter','latex','fontsize',14)
hold off

subplot(2,2,3)
hold on
plot(vTime,100 * vAggregateConsumption_full,'linewidth',1.5,'linestyle','-','color',[.03,.24,.8])
plot(vTime,100 * vAggregateConsumption_reduced,'linewidth',1.5,'linestyle','--','color',[.8,.24,.03])
set(gcf,'color','w')
xlim([vTime(1) vTime(end)])
grid on
title('Consumption','interpreter','latex','fontsize',14)
ylabel('$\%$ deviation from s.s.','interpreter','latex')
xlabel('Quarters','interpreter','latex')
hold off

subplot(2,2,4)
hold on
plot(vTime,100 * vAggregateInvestment_full,'linewidth',1.5,'linestyle','-','color',[.03,.24,.8])
plot(vTime,100 * vAggregateInvestment_reduced,'linewidth',1.5,'linestyle','--','color',[.8,.24,.03])
set(gcf,'color','w')
xlim([vTime(1) vTime(end)])
grid on
title('Investment','interpreter','latex','fontsize',14)
xlabel('Quarters','interpreter','latex')
hold off