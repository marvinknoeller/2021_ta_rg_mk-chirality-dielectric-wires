
N = 5;
for nn=1:N
    NVec(nn^2:nn^2+2*nn)=nn;
end
Q = 2*N*(N+2);
Qd2 = Q/2;
p_big_N = pp(:,:,end);
Var.Q = Q;
Var.NVec = NVec;
epsvec = 1:0.5:50;
parfor kk = 1:length(epsvec)
    eps_rel = epsvec(kk);
    pfinal_sim = pp(:,:,end);
    FF = FarFieldMatrixFunction_Spline_epsN(pfinal_sim,Var,eps_rel);
    [chir_f(kk),smooth_f(kk),cint(kk)] = chiral(FF);
end
f=figure;
f.Position = [2230 753 334 345];
plot(epsvec,chir_f,'--k',epsvec,smooth_f,'-b','LineWidth',3.);
hold on
plot(5*ones(1,10), linspace(0,0.3,10),'-r','LineWidth',2.);
xlim([1,50]);
ylim([0,0.3])
set(gca,'FontSize', 20);
ell = legend('$J_2$','$J_{\rm{HS}}$','Interpreter','Latex',...
    'Fontsize',18,'Location','NorthWest');
ell.Position = [0.25,0.74,0.17 0.17];
grid on
ax = gca;
ax.GridAlpha = 0.9;
ax.YMinorGrid = "off";
ax.MinorGridAlpha = 0.6;
ax.MinorGridLineStyle = "--";
ax.LineWidth = 1.2;
set(gca,'FontSize', 14);
xlabel('$\varepsilon_r$','Interpreter','Latex','Fontsize',20);
set(gca,'GridAlpha', 0.5);
filename = "PaperPlots/eps_meas1";
print(gcf,'-depsc',filename);
