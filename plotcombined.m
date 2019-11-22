function grafica = plotcombined(x,y,tituloFigura,formato,nFig)
% x = k\delta x
% y = M o k*\delta x
% titulo es el titulo de la grafica
% tituloFigura es el nombre para el archivo
% formato si se trata de la amplitud (A) o de la fase (F)
grafica = figure(nFig)
for i = 1:size(y,1)
    plot(x,y(i,:))
    hold on
end
grid on
set(legend,'Interpreter','latex')
xlabel('k\Delta x')
legend('FOU','UCW3-ADER3','UCW3-ADER2','UCW3-ADER1','CD','UPW1-RK2','UCW3-RK2','UCW3-RK3','Location','northeast')
if formato == 'A'
    set(gca, 'YScale', 'log')
    ylabel('\eta/ak')
    xlim([0 2*pi])
    set(gca,'XTickLabel',{'0','\pi/2','\pi','3\pi/2','2\pi'})
    saveas(grafica,strcat(tituloFigura, '_diffusion'),'epsc')
else
    ylabel('k*\Delta x')
    xlim([0 pi])
    set(gca,'XTick',0:pi/4:pi) 
    set(gca,'XTickLabel',{'0','\pi/4','\pi/2','3\pi/4','\pi'})
    saveas(grafica,strcat(tituloFigura, '_dispersion'),'epsc')
end
hold off
end