steps = size(pp,3)-1;
f=figure;
f.Position = [2230 753 334 345];
plot(0:steps-1,chir(1:steps),'--k',0:steps-1,smooth_relax(1:steps),'-b','LineWidth',3);
xlim([0,steps-1]);
ylim([0,0.4])
legend('$J_2$','$J_{\rm{HS}}$','Interpreter','Latex',...
    'Fontsize',18,'Location','NorthWest');

grid on
ax = gca;
ax.GridAlpha = 0.9;
ax.YMinorGrid = "off";
ax.MinorGridAlpha = 0.6;
ax.MinorGridLineStyle = "--";
ax.XAxis.TickValues = 0:20:steps;
ax.LineWidth = 1.2;
set(gca,'FontSize', 14);
xlabel('$\ell$','Interpreter','Latex','Fontsize',20);
set(gca,'GridAlpha', 0.5);
filename = "PaperPlots/iter_meas2";
print(gcf,'-depsc',filename);