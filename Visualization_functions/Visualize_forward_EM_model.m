% Emission Model
% Divya Kumawat, 09/2022
%% Visualizing the forward model

function Visualize_forward_EM_model(wc_soil,tbh_tauomega,tbv_tauomega)

cmap(1,:) = 'r';
cmap(2,:) = 'b';

h = figure;
set(h, 'Position', [337,310,1318,398])
subplot(1,3,1)
plot(wc_soil,tbh_tauomega(:,1,1), 'Color',cmap(1,:),'linewidth',1.5,...
    'Marker', 'o','MarkerSize',6,'MarkerIndices',1:15:100,MarkerFaceColor=cmap(1,:)); hold on;
plot(wc_soil,tbv_tauomega(:,1,1), 'Color',cmap(1,:),'LineStyle','--','linewidth',1.5,...
    'Marker', 'o','MarkerSize',6,'MarkerIndices',1:15:100,MarkerFaceColor=cmap(1,:)); hold on;

plot(wc_soil,tbh_tauomega(:,1,2), 'Color',cmap(2,:),'linewidth',1.5,...
    'Marker', '^','MarkerSize',6,'MarkerIndices',1:15:100,MarkerFaceColor=cmap(2,:)); hold on;
plot(wc_soil,tbv_tauomega(:,1,2), 'Color',cmap(2,:),'LineStyle','--','linewidth',1.5,...
    'Marker', '^','MarkerSize',6,'MarkerIndices',1:15:100,MarkerFaceColor=cmap(2,:)); hold on;
xlabel("$\theta$ [$\rm m^3.m^{-3}$]",'fontsize',18,"Interpreter","latex")
ylabel("$T_b$ [$\rm K$]",'fontsize',18,"Interpreter","latex")
title('VOD = 0','fontsize',18,"Interpreter","latex")
set(gca,'fontsize',15,'LineWidth', 2,'TickLabelInterpreter','latex');
xlim([0.01,0.6])
ylim([110,280])
grid on

subplot(1,3,2)
plot(wc_soil,tbh_tauomega(:,2,1), 'Color',cmap(1,:),'linewidth',1.5,...
    'Marker', 'o','MarkerSize',6,'MarkerIndices',1:15:100,MarkerFaceColor=cmap(1,:)); hold on;
plot(wc_soil,tbv_tauomega(:,2,1), 'Color',cmap(1,:),'LineStyle','--','linewidth',1.5,...
    'Marker', 'o','MarkerSize',6,'MarkerIndices',1:15:100,MarkerFaceColor=cmap(1,:)); hold on;

plot(wc_soil,tbh_tauomega(:,2,2), 'Color',cmap(2,:),'linewidth',1.5,...
    'Marker', '^','MarkerSize',6,'MarkerIndices',1:15:100,MarkerFaceColor=cmap(2,:)); hold on;
plot(wc_soil,tbv_tauomega(:,2,2), 'Color',cmap(2,:),'LineStyle','--','linewidth',1.5,...
    'Marker', '^','MarkerSize',6,'MarkerIndices',1:15:100,MarkerFaceColor=cmap(2,:)); hold on;
xlabel("$\theta$ [$\rm m^3.m^{-3}$]",'fontsize',18,"Interpreter","latex")
title('VOD = 0.25','fontsize',18,"Interpreter","latex")
set(gca,'fontsize',15,'LineWidth', 2,'TickLabelInterpreter','latex');
xlim([0.01,0.6])
ylim([180,270])
grid on

subplot(1,3,3)
plot(wc_soil,tbh_tauomega(:,3,1), 'Color',cmap(1,:),'linewidth',1.5,...
    'Marker', 'o','MarkerSize',6,'MarkerIndices',1:15:100,MarkerFaceColor=cmap(1,:)); hold on;
plot(wc_soil,tbv_tauomega(:,3,1), 'Color',cmap(1,:),'LineStyle','--','linewidth',1.5,...
    'Marker', 'o','MarkerSize',6,'MarkerIndices',1:15:100,MarkerFaceColor=cmap(1,:)); hold on;

plot(wc_soil,tbh_tauomega(:,3,2), 'Color',cmap(2,:),'linewidth',1.5,...
    'Marker', '^','MarkerSize',6,'MarkerIndices',1:15:100,MarkerFaceColor=cmap(2,:)); hold on;
plot(wc_soil,tbv_tauomega(:,3,2), 'Color',cmap(2,:),'LineStyle','--','linewidth',1.5,...
    'Marker', '^','MarkerSize',6,'MarkerIndices',1:15:100,MarkerFaceColor=cmap(2,:)); hold on;
xlabel("$\theta$ [$\rm m^3.m^{-3}$]",'fontsize',18,"Interpreter","latex")
title('VOD = 0.5','fontsize',18,"Interpreter","latex")
grid on
set(gca,'fontsize',15,'LineWidth', 2,'TickLabelInterpreter','latex');
xlim([0.01,0.6])
ylim([220,270])
dummyh11 = plot(nan, nan, 'Linestyle', '-', 'Marker', 'none', 'Color', 'k','linewidth',1.5);
dummyh21 = plot(nan, nan, 'Linestyle', '--', 'Marker', 'none', 'Color', 'k','linewidth',1.5);
hold on;
dummyh1 = scatter(nan, nan,'filled', 'Marker', 'o', 'MarkerFaceColor', cmap(1,:),'linewidth',1.5);
dummyh2 = scatter(nan, nan, 'filled', 'Marker', '^', 'MarkerFaceColor', cmap(2,:),'linewidth',1.5);
legend([dummyh11,dummyh21,dummyh1,dummyh2],{'H-pol','V-pol', '$\rho_s$ = 100 ','$\rho_s$ = 400 '},...
    "Interpreter","latex","Position",[0.11,0.23,0.14183,0.29015],"FontSize",11, 'EdgeColor','none',...
    'Color','none');
end