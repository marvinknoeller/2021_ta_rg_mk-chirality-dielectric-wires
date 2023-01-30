choice = [1,11,31, 51, 71, size(pp,3)];
for kk = 1: length(choice)
    iteration = choice(kk);
    f=figure;
    f.Position = [680 753 334 345];
    f.Color = 'W';
    X_stars = pp(:,:,iteration);
    X = splinepoints(pp(:,:,iteration),11);
    plot3(X(1,:),X(2,:),X(3,:),'-','LineWidth', 1,'Color','b', 'LineWidth', 2.5)
    hold all
    plot3(ones(size(X(1,:)))*4,X(2,:),X(3,:),'-k', 'LineWidth', 2.5)
    plot3(X(1,:),ones(size(X(2,:)))*4.001,X(3,:),'-k', 'LineWidth', 2.5)
    plot3(X(1,:),X(2,:),ones(size(X(3,:)))*(-4),'-k', 'LineWidth', 2.5)
    hold off
    axis equal
    axis([-4.,4.,-4.001,4.001,-4,4])
    set(gca,'fontsize',14)
    grid on
    set(gca,'GridAlpha', 0.2);
    set(gca,'LineWidth',2.,'TickLength',[0.025 0.04]);
    ax = gca;
    ax.XAxis.TickValues = [-4,-2,0,2,4];
    ax.ZAxis.TickValues = [-4,-2,0,2,4];
    ax.XAxis.TickLabels(5) = {''};
    hold off
    rem = "";
    if iteration == size(pp,3)
        form = 2;
        for2 = 1;
        iterationnum = iteration-1;
        rem = " (final result)";
    elseif (iteration == 1)
        form = '%10.0e';
        for2 = form;
        iterationnum = iteration;
        rem = " (initial curve)";
    elseif (iteration == choice(2))
        form = '%10.0e';
        iterationnum = iteration;
        for2 = form;
    else
        iterationnum = iteration;
        form = 2;
        for2 = 1;
    end
    title(strcat("$\ell=$",num2str(iteration-1),rem),'Interpreter','Latex','FontSize',18)
    text(3.11,-7.31,-1.21,{...
        strcat('$$J_{\rm{HS}}=$$',num2str(smooth_relax(iterationnum),for2))}...
        ,'BackGroundColor','w',...
            'FontSize',18,'Interpreter','Latex',...
            'Edgecolor','k',...
            'LineStyle','-',...
            'LineWidth',1,...
            'Clipping','off',...
            'HorizontalAlignment', 'right',...
            'VerticalAlignment', 'top')
        filename=strcat("PaperPlots/1Iterates",num2str(kk));
        print(gcf,'-depsc',filename);
end