workers10 = hdf5read('../../test//reference//timingAndWorkers.h5', 'times10');
workers30 = hdf5read('../../test//reference//timingAndWorkers30.h5', 'times30');

%figure('Position',[20,20,900,600],'Name',...
%    'Log Likelihood Increments: TPF vs KF','Color','w')

plot(workers10,'LineStyle','-','Color','g','LineWidth',2.5) % from Kalman filter
hold on
plot(workers30,'LineStyle','-','Color','b','LineWidth',2.5) % from Matlab filter
%hold on 
xlabel('Time steps')
ylabel('Time per mutation(s)')
legend('10 Workers','30 Workers', 'Location', 'southwest')
%axis([USquarter(1) USquarter(end) -25 max(max(lik_TPF), max(lik_TPF))+1])
%set(gca,'FontSize',20)