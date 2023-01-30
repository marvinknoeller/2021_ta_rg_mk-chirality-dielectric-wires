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
    plot3(ones(size(X(1,:)))*7,X(2,:),X(3,:),'-k', 'LineWidth', 2.5)
    plot3(X(1,:),ones(size(X(2,:)))*7.001,X(3,:),'-k', 'LineWidth', 2.5)
    plot3(X(1,:),X(2,:),ones(size(X(3,:)))*(-7),'-k', 'LineWidth', 2.5)
    hold off
    axis equal
    axis([-7.,7.,-7.001,7.001,-7,7])
    set(gca,'fontsize',14)
    grid on
    set(gca,'GridAlpha', 0.2);
    set(gca,'LineWidth',2.,'TickLength',[0.025 0.04]);
    ax = gca;
    ax.XAxis.TickValues = [-7,-5,-3,-1,1,3,5,7];
    ax.YAxis.TickValues = [-7,-5,-3,-1,1,3,5,7];
    ax.ZAxis.TickValues = [-7,-5,-3,-1,1,3,5,7];
    ax.XAxis.TickLabels(7) = {''};
    ax.XAxis.TickLabels(8) = {''};
    hold off
    rem = "";
    if iteration == size(pp,3)
        form = 2;
        for2 = 2;
        iterationnum = iteration-1;
        rem = " (final result)";
        add_zero = '';
    elseif (iteration == 1)
        form = '%10.0e';
        for2 = form;
        iterationnum = iteration;
        rem = " (initial curve)";
        add_zero = '';
    elseif (iteration == choice(2)) 
        iterationnum = iteration;
        form = 2;
        for2 = 1;
        add_zero = '0';
    else
        iterationnum = iteration;
        form = 2;
        for2 = 2;
        add_zero = '';
    end
    title(strcat("$\ell=$",num2str(iteration-1),rem),'Interpreter','Latex','FontSize',18)

text(7.31,-9.99,-4.66,{...
        strcat('$$J_{\rm{HS}}=$$',num2str(smooth_relax(iterationnum),for2))}...
        ,'BackGroundColor','w',...
            'FontSize',18,'Interpreter','Latex',...
            'Edgecolor','k',...
            'LineStyle','-',...
            'LineWidth',1,...
            'Clipping','off',...
            'HorizontalAlignment', 'right',...
            'VerticalAlignment', 'top')
        filename=strcat("PaperPlots/2Iterates",num2str(kk));
        print(gcf,'-depsc',filename);
end