% Emission Model
% Divya Kumawat, 09/2022
%% Plots density scatter plot for synthetic retrievals

function Scatter_plot(combo,xopt)

error=combo(:,1:2)-xopt;

[R2,NSE,RMSE,Bias]=percriteria(combo(:,1),xopt(:,1));
Performancetau=[std(error(:,1)),R2,NSE,RMSE,Bias];

[R2,NSE,RMSE,Bias]=percriteria(combo(:,2),xopt(:,2));
PerformanceRos=[std(error(:,2)),R2,NSE,RMSE,Bias];

figure;
subplot(1,2,1);
scatter_kde(combo(:,1),xopt(:,1),'filled', 'MarkerSize', 12);
load('colormap_dscatter4.mat');
colormap(colormapdscatter4);
grid on
box on
hline = refline(1,0);
hline.Color = 'r';
hline.LineWidth = 1;
box on
ax = gca;
ax.LineWidth = 1;
set(gca,'FontSize',14);
xlabel('$\tau$','FontSize',20,'FontWeight','bold','Interpreter','Latex')
ylabel('$$\hat{\tau}$$','FontSize',20,'FontWeight','bold','Interpreter','Latex')
yylim=get(gca,'ylim');xxlim=get(gca,'xlim');
text(0.49,0.11,0,{['$r^2$ = ',num2str(Performancetau(1,2),'%.3f')]
    ['Bias = ',num2str(Performancetau(1,5),'%.3f')]
    ['RMSE = ',num2str(Performancetau(1,4),'%.3f')]},'FontSize',14,'HorizontalAlignment','right','Interpreter','Latex');
set(gca,'TickLabelInterpreter','latex');
axis square;
xlim([0,0.5]);
ylim([0,0.5]);


subplot(1,2,2)
scatter_kde(combo(:,2),xopt(:,2),'filled', 'MarkerSize', 12);
load('colormap_dscatter4.mat');
colormap(colormapdscatter4);
grid on
box on
hline = refline(1,0);
hline.Color = 'r';
hline.LineWidth = 1;
box on
ax = gca;
ax.LineWidth = 1;
set(gca,'FontSize',14)
xlabel('$\theta $','FontSize',14,'Interpreter','Latex')
ylabel('$$\hat{\theta}$$','FontSize',14,'Interpreter','Latex')
yylim=get(gca,'ylim');xxlim=get(gca,'xlim');
text(0.6,0.13,0,{['$r^2$ = ',num2str(PerformanceRos(1,2),'%.3f')]
    ['Bias = ',num2str(PerformanceRos(1,5),'%.3f')]
    ['RMSE = ',num2str(PerformanceRos(1,4),'%.3f')]},'FontSize',14,'HorizontalAlignment','right','Interpreter','Latex');
set(gca,'TickLabelInterpreter','latex')
axis square
xlim([0,0.6])
ylim([0,0.6])
end