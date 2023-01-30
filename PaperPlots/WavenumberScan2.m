N = 14;
for nn=1:N
    NVec(nn^2:nn^2+2*nn)=nn;
end
Q = 2*N*(N+2);
Qd2 = Q/2;
p_big_N = pp(:,:,end);
Var.Q = Q;
Var.NVec = NVec;
kappaVec = linspace(0.2,2,50);
parfor kk = 1:length(kappaVec)
    kappa = kappaVec(kk);
    pfinal_sim = pp(:,:,end);
    FF = FarFieldMatrixFunction_Spline_kappaN(pfinal_sim,Var,kappa);
    [chir_f(kk),smooth_f(kk),cint(kk)] = chiral(FF);
end
f=figure;
f.Position = [2230 753 334 345];
plot(kappaVec,chir_f,'--k',kappaVec,smooth_f,'-b','LineWidth',3.);
hold on
plot(1*ones(1,10), linspace(0,0.4,10),'-r','LineWidth',2.);
xlim([0.2,2]);
ylim([0,0.4])
set(gca,'FontSize', 20);
legend('$J_2$','$J_{\rm{HS}}$','Interpreter','Latex',...
    'Fontsize',18,'Location','NorthWest');
grid on
ax = gca;
ax.GridAlpha = 0.9;
ax.YMinorGrid = "off";
ax.MinorGridAlpha = 0.6;
ax.MinorGridLineStyle = "--";
ax.LineWidth = 1.2;
set(gca,'FontSize', 14);
xlabel('$k$','Interpreter','Latex','Fontsize',20);
set(gca,'GridAlpha', 0.5);
filename = "PaperPlots/k_meas2";
print(gcf,'-depsc',filename);
