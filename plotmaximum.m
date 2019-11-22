function grafica = plotmaximum(x,y,tituloFigura,nFig)
% x = k\delta x
% y = M o k*\delta x
% titulo es el titulo de la grafica
% tituloFigura es el nombre para el archivo
% formato si se trata de la amplitud (A) o de la fase (F)
grafica = figure(nFig)
color = distinguishable_colors(length(x));

for i = 1:size(y,1)
    plot(x,y(i,:),'LineWidth',2,'color',color(i,:))
    hold on
end
grid on
set(legend,'Interpreter','latex')
xlabel('\boldmath$$\sigma$$','Interpreter','latex');
legend('FOU','CD','UCW3-ADER1','UCW3-ADER2','UCW3-ADER3','UPW1-RK2','UCW3-RK2','UCW3-RK3','Location','northwest')
ylabel('\boldmath$$|M|_{max}$$','Interpreter','latex');
xlim([0 1.5])
xticks([0 0.25 0.5 0.75 1.0 1.25 1.5]);
ylim([1 3])
saveas(grafica,strcat(tituloFigura, '_maximum'),'epsc')
hold off
end