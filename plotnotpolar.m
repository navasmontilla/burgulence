function grafica = plotnotpolar(x,y,titulo,tituloFigura,formato,nFig)
% x = k\delta x
% y = M o k*\delta x
% titulo es el titulo de la grafica
% tituloFigura es el nombre para el archivo
% formato si se trata de la amplitud (A) o de la fase (F)
grafica = figure(nFig)
hold on
title(titulo)

for i = 1:floor(size(y,1)/2)
    plot(x,y(i,:),'LineWidth',1.5)
end
for i = floor(size(y,1)/2)+1:size(y,1)
    plot(x,y(i,:),'--','LineWidth',1.5)
end
grid on


xlabel('\boldmath$$k\Delta x$$','Interpreter','latex');
switch formato
    case 'A'
        ylabel('\boldmath$$|M|$$','Interpreter','latex');
        xlim([0 2*pi]);
        % ylim([0 1]);
        set(gca,'XTick',0:pi/2:2*pi);
        set(gca,'XTickLabel',{'0','\pi/2','\pi','3\pi/2','2\pi'});
        saveas(grafica,strcat(tituloFigura, '_diffusion'),'epsc');
        set(legend,'Interpreter','latex')
        if size(y,1)>1

            lgd = legend('0.0','0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','1.0','Location','southwest')
            lgd.Title.String = 'CFL';
        end
        lgd = legend('0.001','0.100','0.200','0.300','0.400','0.500','0.600','0.700','0.800','0.900','1.000','Location','southwest');
        lgd.Title.String = 'CFL';
    case 'F' 
        plot(x,x,':','Color','black','LineWidth',1)
        ylabel('\boldmath$$\tilde{k}\Delta x$$','Interpreter','latex');
        xlim([0 pi])
        ylim([-pi pi]);
        set(gca,'XTick',0:pi/4:pi);
        set(gca,'XTickLabel',{'0','\pi/4','\pi/2','3\pi/4','\pi'});
        set(gca,'YTick',-pi:pi/4:pi);
        set(gca,'YTickLabel',{'-\pi','-3\pi/4','-\pi/2','-\pi/4','0','\pi/4','\pi/2','3\pi/4','\pi'});
        saveas(grafica,strcat(tituloFigura, '_dispersion'),'epsc')
    case 'Z'
        ylabel('\boldmath$$|M|$$','Interpreter','latex');
        xlim([0 pi/2])
        ylim([0.7 1.0]);
        set(gca,'XTick',0:pi/8:pi/2) 
        set(gca,'XTickLabel',{'0','\pi/8','\pi/4','3\pi/8','\pi/2'})
        saveas(grafica,strcat(tituloFigura, '_diffusion_zoom'),'epsc')
        lgd = legend('0.001','0.100','0.200','0.300','0.400','0.500','0.600','0.700','0.800','0.900','1.000','Location','southwest');
        lgd.Title.String = 'CFL';
end
hold off
end