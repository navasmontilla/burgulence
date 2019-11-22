function grafica = plotpolar(x,y,titulo,tituloFigura,nFig)
% x = k\delta x
% y = M o k*\delta x
% titulo es el titulo de la grafica
% tituloFigura es el nombre para el archivo
grafica = figure(nFig)
polaraxes
title(titulo)
for i = 1:floor(size(y,1)/2)
    polarplot(x,y(i,:),'LineWidth',1.5)
    hold on
end
for i = floor(size(y,1)/2)+1:size(y,1)
    polarplot(x,y(i,:),'--','LineWidth',1.5)
    hold on
end
grid on
rlim([0 1]);
rticks(0:0.5:2.5)
lgd = legend('0.0','0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','1.0','Location','southeast')
lgd.Title.String = 'CFL';
hold off
saveas(grafica,strcat(tituloFigura, '_diffusion_polar'),'epsc')
end