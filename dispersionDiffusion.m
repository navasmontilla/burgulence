%% Analisis de estabilidad
clear all;
close all;
%% Variables globales
nu          = [0.001 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0];
betaVector  = 0:0.001:2*pi;
CFLVector   = 0:0.001:1.5;
limite      = 100.0*ones(8,1);
%% First Order Upwind
for i = 1:11
    M(i,:)       = sqrt((1-nu(i)+nu(i).*cos(betaVector)).^2+(nu(i).*sin(betaVector)).^2);
    theta(i,:)   = atan((nu(i).*sin(betaVector))./(1+nu(i).*cos(betaVector)-nu(i)));
    betaMod(i,:) = theta(i,:)./nu(i) ;
end
% Semidiscreto
semi_betaMod(1,:)     = sin(betaVector);
semi_M(1,:)           = (1-cos(betaVector))./betaVector;
limite(1)             = 100.0;
for i = 1:length(CFLVector)
    Mmax(1,i)      = max(sqrt((1-CFLVector(i)+CFLVector(i).*cos(betaVector)).^2+(CFLVector(i).*sin(betaVector)).^2));
    if (Mmax(1,i) > 1.001) && (CFLVector(i) < limite(1,1))
        limite(1,1) = CFLVector(i);
    end
end
% Representacion
titulo           = 'FOU';
tituloFigura     = 'FOU';
plotnotpolar(betaVector,M,titulo,tituloFigura,'A',1)
plotnotpolar(betaVector,betaMod,titulo,tituloFigura,'F',2)
plotpolar(betaVector,M,titulo,tituloFigura,3)
%% Central Difference
for i = 1:11
   M(i,:)       = 1-nu(i).*sin(betaVector);
   theta(i,:)   = 2*pi./nu(i)
   betaMod(i,:) = theta(i,:)./nu(i);
end
semi_betaMod(2,:)     = sin(betaVector);
semi_M(2,:)           = zeros(size(semi_M,1));
for i = 1:length(CFLVector)
    Mmax(2,i)      = max(1-CFLVector(i).*sin(betaVector));
    if (Mmax(2,i) > 1.001) && (CFLVector(i) < limite(2,1))
        limite(2,1) = CFLVector(i);
    end
end
% Representacion
titulo = 'CD';
tituloFigura = 'CD';
plotnotpolar(betaVector,M,titulo,tituloFigura,'A',4)
plotnotpolar(betaVector,betaMod,titulo,tituloFigura,'F',5)
plotpolar(betaVector,M,titulo,tituloFigura,6)
%% Third Order Upwind - ADER 1
for i = 1:11
    M(i,:)       = sqrt((1-1/2*nu(i)+2./3*nu(i).*cos(betaVector)-1./6*nu(i).*cos(2.*betaVector)).^2+(4./3*nu(i).*sin(betaVector)-1./6*nu(i).*sin(2.*betaVector)).^2);
    theta(i,:)   = atan((4/3*nu(i).*sin(betaVector)-1/6*nu(i).*sin(2.*betaVector))./(1-1/2*nu(i)+2/3*nu(i).*cos(betaVector)-1/6*nu(i).*cos(2.*betaVector)));
    betaMod(i,:) = theta(i,:)./nu(i) ;
end
semi_betaMod(3,:)     = 4./3.*sin(betaVector)-1./6.*sin(2.*betaVector);
semi_M(3,:)           = (-1./2-1/6.*cos(2.*betaVector)-2/3.*cos(betaVector))./betaVector;
for i = 1:length(CFLVector)
    Mmax(3,i)      = max(sqrt((1-1/2*CFLVector(i)+2./3*CFLVector(i).*cos(betaVector)-1./6*CFLVector(i).*cos(2.*betaVector)).^2+(4./3*CFLVector(i).*sin(betaVector)-1./6*CFLVector(i).*sin(2.*betaVector)).^2));
    if (Mmax(3,i) > 1.001) && (CFLVector(i) < limite(3,1))
        limite(3,1) = CFLVector(i);
    end
end
% Representacion
titulo = 'UCW3-ADER1';
tituloFigura = 'UCW3-ADER1';
plotnotpolar(betaVector,M,titulo,tituloFigura,'A',7)
plotnotpolar(betaVector,betaMod,titulo,tituloFigura,'F',8)
plotpolar(betaVector,M,titulo,tituloFigura,9)
%% Third Order Upwind - ADER 2
for i = 1:11
    M(i,:)       = sqrt((1-1/2*nu(i)+2./3*nu(i).*cos(betaVector)-1./6*nu(i).*cos(2.*betaVector)+nu(i)^2.*(cos(betaVector)-1)).^2+(4./3*nu(i).*sin(betaVector)-1./6*nu(i).*sin(2.*betaVector)).^2);
    theta(i,:)   = atan((4/3*nu(i).*sin(betaVector)-1/6*nu(i).*sin(2.*betaVector))./(1-1/2*nu(i)+2/3*nu(i).*cos(betaVector)-1/6*nu(i).*cos(2.*betaVector)+nu(i)^2.*(cos(betaVector)-1)));
    betaMod(i,:) = theta(i,:)./nu(i);
end
% DUDAS!
semi_betaMod(4,:)     = 4./3.*sin(betaVector)-1./6.*sin(2.*betaVector);
semi_M(4,:)           = (-1./2-1/6.*cos(2.*betaVector)-2/3.*cos(betaVector))./betaVector;
for i = 1:length(CFLVector)
    Mmax(4,i)      = max(sqrt((1-1/2*CFLVector(i)+2./3*CFLVector(i).*cos(betaVector)-1./6*CFLVector(i).*cos(2.*betaVector)+CFLVector(i)^2.*(cos(betaVector)-1)).^2+(4./3*CFLVector(i).*sin(betaVector)-1./6*CFLVector(i).*sin(2.*betaVector)).^2));
    if (Mmax(4,i) > 1.001) && (CFLVector(i) < limite(4,1))
        limite(4,1) = CFLVector(i);
    end
end
% Representacion
titulo = 'UCW3-ADER2';
tituloFigura = 'UCW3-ADER2';
plotnotpolar(betaVector,M,titulo,tituloFigura,'A',10)
plotnotpolar(betaVector,betaMod,titulo,tituloFigura,'F',11)
plotpolar(betaVector,M,titulo,tituloFigura,12)
%% Third Order Upwind - ADER 3
for i = 1:11
    M(i,:)       = sqrt((1-nu(i)/2*(1+2*nu(i)-nu(i)^2)+nu(i).*cos(betaVector).*(2/3+nu(i)-2/3*nu(i)^2)+nu(i)/6.*cos(2.*betaVector).*(nu(i)^2-1)).^2+(1/6*nu(i).*sin(2.*betaVector).*(1-nu(i)^2)-nu(i)/3.*sin(betaVector).*(4-nu(i)^2)).^2);
    theta(i,:)   = atan((-1/6*nu(i).*sin(2.*betaVector).*(1-nu(i)^2)+nu(i)/3.*sin(betaVector).*(4-nu(i)^2))./(1-nu(i)/2*(1+2*nu(i)-nu(i)^2)+nu(i).*cos(betaVector).*(2/3+nu(i)-2/3*nu(i)^2)+nu(i)/6.*cos(2.*betaVector).*(nu(i)^2-1)));
    betaMod(i,:) = theta(i,:)./nu(i);
end
% DUDAS!
semi_betaMod(5,:)     = 4./3.*sin(betaVector)-1./6.*sin(2.*betaVector);
semi_M(5,:)           = (-1./2-1/6.*cos(2.*betaVector)+2/3.*cos(betaVector))./betaVector;
for i = 1:length(CFLVector)
    Mmax(5,i)      = max(sqrt((1-CFLVector(i)/2*(1+2*CFLVector(i)-CFLVector(i)^2)+CFLVector(i).*cos(betaVector).*(2/3+CFLVector(i)-2/3*CFLVector(i)^2)+CFLVector(i)/6.*cos(2.*betaVector).*(CFLVector(i)^2-1)).^2+(1/6*CFLVector(i).*sin(2.*betaVector).*(1-CFLVector(i)^2)-CFLVector(i)/3.*sin(betaVector).*(4-CFLVector(i)^2)).^2));
    if (Mmax(5,i) > 1.001) && (CFLVector(i) < limite(5,1))
        limite(5,1) = CFLVector(i);
    end
end
% Representacion
titulo = 'UCW3-ADER3';
tituloFigura = 'UCW3-ADER3';
plotnotpolar(betaVector,M,titulo,tituloFigura,'A',13)
plotnotpolar(betaVector,betaMod,titulo,tituloFigura,'F',14)
plotpolar(betaVector,M,titulo,tituloFigura,15)
% Representacion del semidiscreto
        figure113 = figure(113);
        plot(betaVector,semi_betaMod(1,:),'LineWidth',2,'Color','blue');
        hold on;
        plot(betaVector,semi_betaMod(5,:),'LineWidth',2,'Color','red');
        grid on;
        line([0 2*pi],[0 2*pi],'LineStyle',':','Color','black','LineWidth',0.5);
        legend('FOU','UWC3','Ideal');
        % legend('FOU','Ideal');
        xlim([0 pi]);
        ylim([0 pi/2]);
        xticks([0 pi/8 2*pi/8 3*pi/8 4*pi/8 5*pi/8 6*pi/8 7*pi/8 8*pi/8]);
        xticklabels({'0' '\pi/8' '\pi/4' '3\pi/8' '\pi/2' '5\pi/8' '3\pi/4' '7\pi/8' '\pi'});
        yticks([0:pi/8:pi/2]);
        yticklabels({'0' '\pi/8' '\pi/4' '3\pi/8' '\pi/2'});
        xlabel('\boldmath$$k\Delta x$$','Interpreter','latex');
        ylabel('\boldmath$$\tilde{k}\Delta x$$','Interpreter','latex');
        hold off;
        saveas(figure113,sprintf('dispersion_semidiscret_%s',titulo),'epsc');
        
        figure114 = figure(114);
        plot(betaVector,-semi_M(1,:),'LineWidth',2,'Color','blue');
        hold on;
        plot(betaVector,semi_M(5,:),'LineWidth',2,'Color','red');
        grid on;
        ylabel('\boldmath$$\psi$$','Interpreter','latex');
        xlabel('\boldmath$$k\Delta x$$','Interpreter','latex');
        legend('FOU','UWC3');
        % legend('FOU');
        xlim([0 pi]);
        ylim([-0.8 0]);
        xticks([0 pi/8 2*pi/8 3*pi/8 4*pi/8 5*pi/8 6*pi/8 7*pi/8 8*pi/8]);
        yticks([-0.8:0.2:0]);
        xticklabels({'0' '\pi/8' '\pi/4' '3\pi/8' '\pi/2' '5\pi/8' '3\pi/4' '7\pi/8' '\pi'});
        hold off;
        saveas(figure114,sprintf('diffusion_semidiscret_%s',titulo),'epsc');
%% UPW1 - RK2
for i = 1:11
    M(i,:)       = sqrt((1-nu(i)+nu(i)*nu(i)/2+nu(i).*(1-nu(i)).*cos(betaVector)+nu(i).*nu(i).*cos(2.*betaVector)./2).^2+(nu(i).*(1-nu(i)).*sin(betaVector)+nu(i).*nu(i).*sin(2.*betaVector)./2).^2);
    theta(i,:)   = atan((nu(i).*(1-nu(i)).*sin(betaVector)+nu(i).*nu(i).*sin(2.*betaVector)./2)./(1-nu(i)+nu(i)*nu(i)/2+nu(i).*(1-nu(i)).*cos(betaVector)+nu(i).*nu(i).*cos(2.*betaVector)./2));
    betaMod(i,:) = theta(i,:)./nu(i);
end
semi_betaMod(6,:)     = sin(betaVector);
semi_M(6,:)           = (1-cos(betaVector))./betaVector;
for i = 1:length(CFLVector)
    Mmax(6,i)      = max(sqrt((1-CFLVector(i)+CFLVector(i)*CFLVector(i)/2+CFLVector(i).*(1-CFLVector(i)).*cos(betaVector)+CFLVector(i).*CFLVector(i).*cos(2.*betaVector)./2).^2+(CFLVector(i).*(1-CFLVector(i)).*sin(betaVector)+CFLVector(i).*CFLVector(i).*sin(2.*betaVector)./2).^2));
    if (Mmax(6,i) > 1.001) && (CFLVector(i) < limite(6,1))
        limite(6,1) = CFLVector(i);
    end
end
% Representacion
titulo = 'UW1-RK2';
tituloFigura = 'UW1-RK2';
plotnotpolar(betaVector,M,titulo,tituloFigura,'A',16)
plotnotpolar(betaVector,betaMod,titulo,tituloFigura,'F',17)
plotpolar(betaVector,M,titulo,tituloFigura,18)
%% UCW3 - RK2
% real        = (1-nu(i)/2-5*nu(i)^2/24) + (2*nu(i)/3-5*nu(i)^2/18).*cos(betaVector) + (-nu(i)/6+23*nu(i)^2/36).*cos(2.*betaVector) - (nu(i)^2/6).*cos(3.*betaVector) + (nu(i)^2/72).*cos(4.*betaVector)
% imaginaria  = (-4*nu(i)/3+11*nu(i)^2/18).*sin(betaVector) + (nu(i)/6-19*nu(i)^2/36).*sin(2.*betaVector) + (nu(i)^2/6).*sin(3.*betaVector) - (nu(i)^2/72).*sin(4.*betaVector)
for i = 1:11
    M(i,:)       = sqrt(((1-nu(i)/2-5*nu(i)^2/24) + (2*nu(i)/3-5*nu(i)^2/18).*cos(betaVector) + (-nu(i)/6+23*nu(i)^2/36).*cos(2.*betaVector) - (nu(i)^2/6).*cos(3.*betaVector) + (nu(i)^2/72).*cos(4.*betaVector)).^2 + ((-4*nu(i)/3+11*nu(i)^2/18).*sin(betaVector) + (nu(i)/6-19*nu(i)^2/36).*sin(2.*betaVector) + (nu(i)^2/6).*sin(3.*betaVector) - (nu(i)^2/72).*sin(4.*betaVector)).^2);
    theta(i,:)   = atan(-((-4*nu(i)/3+11*nu(i)^2/18).*sin(betaVector) + (nu(i)/6-19*nu(i)^2/36).*sin(2.*betaVector) + (nu(i)^2/6).*sin(3.*betaVector) - (nu(i)^2/72).*sin(4.*betaVector))./((1-nu(i)/2-5*nu(i)^2/24) + (2*nu(i)/3-5*nu(i)^2/18).*cos(betaVector) + (-nu(i)/6+23*nu(i)^2/36).*cos(2.*betaVector) - (nu(i)^2/6).*cos(3.*betaVector) + (nu(i)^2/72).*cos(4.*betaVector)));
    betaMod(i,:) = theta(i,:)./nu(i);
end
semi_betaMod(7,:)     = 4./3.*sin(betaVector)-1./6.*sin(2.*betaVector);
semi_M(7,:)           = (-1./2-1/6.*cos(2.*betaVector)-2/3.*cos(betaVector))./betaVector;
for i = 1:length(CFLVector)
    Mmax(7,i)      = max(sqrt(((1-CFLVector(i)/2-5*CFLVector(i)^2/24) + (2*CFLVector(i)/3-5*CFLVector(i)^2/18).*cos(betaVector) + (-CFLVector(i)/6+23*CFLVector(i)^2/36).*cos(2.*betaVector) - (CFLVector(i)^2/6).*cos(3.*betaVector) + (CFLVector(i)^2/72).*cos(4.*betaVector)).^2 + ((-4*CFLVector(i)/3+11*CFLVector(i)^2/18).*sin(betaVector) + (CFLVector(i)/6-19*CFLVector(i)^2/36).*sin(2.*betaVector) + (CFLVector(i)^2/6).*sin(3.*betaVector) - (CFLVector(i)^2/72).*sin(4.*betaVector)).^2));
    if (Mmax(7,i) > 1.001) && (CFLVector(i) < limite(7,1))
        limite(7,1) = CFLVector(i);
    end
end
% Representacion
titulo = 'UWC3-RK2';
tituloFigura = 'UWC3-RK2';
plotnotpolar(betaVector,M,titulo,tituloFigura,'A',19)
plotnotpolar(betaVector,betaMod,titulo,tituloFigura,'F',20)
plotpolar(betaVector,M,titulo,tituloFigura,21)
%% UCW3 - RK3
% real        = (1-nu(i)/2-5*nu(i)^2/24+59*nu(i)^3/432) + (2*nu(i)/3-5*nu(i)^2/18-4*nu(i)^3/72).*cos(betaVector) + (-nu(i)/6+23*nu(i)^2/36-35*nu(i)^3/144).*cos(2.*betaVector) + (-nu(i)^2/6+155*nu(i)^3/648).*cos(3.*betaVector) + (nu(i)^2/72-13*nu(i)^3/144).*cos(4.*betaVector) + (nu(i)^3/72).*cos(5.*betaVector) - (nu(i)^3/1296).*cos(6.*betaVector)
% imaginaria  = (-4*nu(i)/3+11*nu(i)^2/18+6*nu(i)^3/72).*sin(betaVector) + (nu(i)/6-19*nu(i)^2/36+27*nu(i)^3/144).*sin(2.*betaVector) + (nu(i)^2/6-163*nu(i)^3/648).*sin(3.*betaVector) + (-nu(i)^2/72+13*nu(i)^3/144).*sin(4.*betaVector) - (nu(i)^3/72).*sin(5.*betaVector) + (nu(i)^3/1296).*sin(6.*betaVector)
for i = 1:11
    M(i,:)       = sqrt(((1-nu(i)/2-5*nu(i)^2/24+59*nu(i)^3/432) + (2*nu(i)/3-5*nu(i)^2/18-4*nu(i)^3/72).*cos(betaVector) + (-nu(i)/6+23*nu(i)^2/36-35*nu(i)^3/144).*cos(2.*betaVector) + (-nu(i)^2/6+155*nu(i)^3/648).*cos(3.*betaVector) + (nu(i)^2/72-13*nu(i)^3/144).*cos(4.*betaVector) + (nu(i)^3/72).*cos(5.*betaVector) - (nu(i)^3/1296).*cos(6.*betaVector)).^2+((-4*nu(i)/3+11*nu(i)^2/18+6*nu(i)^3/72).*sin(betaVector) + (nu(i)/6-19*nu(i)^2/36+27*nu(i)^3/144).*sin(2.*betaVector) + (nu(i)^2/6-163*nu(i)^3/648).*sin(3.*betaVector) + (-nu(i)^2/72+13*nu(i)^3/144).*sin(4.*betaVector) - (nu(i)^3/72).*sin(5.*betaVector) + (nu(i)^3/1296).*sin(6.*betaVector)).^2)
    theta(i,:)   = atan(-((-4*nu(i)/3+11*nu(i)^2/18+6*nu(i)^3/72).*sin(betaVector) + (nu(i)/6-19*nu(i)^2/36+27*nu(i)^3/144).*sin(2.*betaVector) + (nu(i)^2/6-163*nu(i)^3/648).*sin(3.*betaVector) + (-nu(i)^2/72+13*nu(i)^3/144).*sin(4.*betaVector) - (nu(i)^3/72).*sin(5.*betaVector) + (nu(i)^3/1296).*sin(6.*betaVector))./((1-nu(i)/2-5*nu(i)^2/24+59*nu(i)^3/432) + (2*nu(i)/3-5*nu(i)^2/18-4*nu(i)^3/72).*cos(betaVector) + (-nu(i)/6+23*nu(i)^2/36-35*nu(i)^3/144).*cos(2.*betaVector) + (-nu(i)^2/6+155*nu(i)^3/648).*cos(3.*betaVector) + (nu(i)^2/72-13*nu(i)^3/144).*cos(4.*betaVector) + (nu(i)^3/72).*cos(5.*betaVector) - (nu(i)^3/1296).*cos(6.*betaVector)))
    betaMod(i,:) = theta(i,:)./nu(i);
end
semi_betaMod(8,:)     = 4./3.*sin(betaVector)-1./6.*sin(2.*betaVector);
semi_M(8,:)           = (-1./2-1/6.*cos(2.*betaVector)-2/3.*cos(betaVector))./betaVector;
for i = 1:length(CFLVector)
    Mmax(8,i)      = max(sqrt(((1-CFLVector(i)/2-5*CFLVector(i)^2/24+59*CFLVector(i)^3/432) + (2*CFLVector(i)/3-5*CFLVector(i)^2/18-4*CFLVector(i)^3/72).*cos(betaVector) + (-CFLVector(i)/6+23*CFLVector(i)^2/36-35*CFLVector(i)^3/144).*cos(2.*betaVector) + (-CFLVector(i)^2/6+155*CFLVector(i)^3/648).*cos(3.*betaVector) + (CFLVector(i)^2/72-13*CFLVector(i)^3/144).*cos(4.*betaVector) + (CFLVector(i)^3/72).*cos(5.*betaVector) - (CFLVector(i)^3/1296).*cos(6.*betaVector)).^2+((-4*CFLVector(i)/3+11*CFLVector(i)^2/18+6*CFLVector(i)^3/72).*sin(betaVector) + (CFLVector(i)/6-19*CFLVector(i)^2/36+27*CFLVector(i)^3/144).*sin(2.*betaVector) + (CFLVector(i)^2/6-163*CFLVector(i)^3/648).*sin(3.*betaVector) + (-CFLVector(i)^2/72+13*CFLVector(i)^3/144).*sin(4.*betaVector) - (CFLVector(i)^3/72).*sin(5.*betaVector) + (CFLVector(i)^3/1296).*sin(6.*betaVector)).^2));
    if (Mmax(8,i) > 1.001) && (CFLVector(i) < limite(8,1))
        limite(8,1) = CFLVector(i);
    end
end
% Representacion
titulo = 'UWC3-RK3';
tituloFigura = 'UWC3-RK3';
plotnotpolar(betaVector,M,titulo,tituloFigura,'A',22)
plotnotpolar(betaVector,betaMod,titulo,tituloFigura,'F',23)
plotpolar(betaVector,M,titulo,tituloFigura,24)
%% Representacion graficas conjuntas
plotcombined(betaVector,semi_M,'semidiscret','A',25)
plotcombined(betaVector,semi_betaMod,'semidiscret','F',26)
plotmaximum(CFLVector,Mmax,'CFL',27)
%% Estudio de cómo afecta el paso temporal - FOU
for i = 1:11
    K(i)         = floor(1/nu(i));
    nuMod(i)     = (1-K(i)*nu(i))/1; 
    M(i,:)       = sqrt((1-nu(i)+nu(i).*cos(betaVector)).^2+(nu(i).*sin(betaVector)).^2);
    MMod(i,:)    = sqrt((1-nuMod(i)+nuMod(i).*cos(betaVector)).^2+(nuMod(i).*sin(betaVector)).^2);
    MTotal(i,:)  = M(i,:).^K(i).*MMod(i,:);
end
titulo = 'FOU';
tituloFigura = 'FOU_DT';
plotnotpolar(betaVector,MTotal,titulo,tituloFigura,'A',28);
plotnotpolar(betaVector,MTotal,titulo,tituloFigura,'Z',29);
% Representacion del número de pasos temporales necesarios
grafica = figure(30)
bar(nu,log10(K));
hold on;
xlabel('CFL');
ylabel('log_{10}(K)');
set(gca,'XTick',0:0.1:1);
grid on;
for i1=1:numel(K)
    text(nu(i1),log10(K(i1)),num2str(K(i1),'%d'),...
               'HorizontalAlignment','center',...
               'VerticalAlignment','bottom')
end
saveas(grafica,strcat(tituloFigura, '_pasosTemporales'),'epsc')
hold off;
%% Estudio de cómo afecta el paso temporal - UCW3-ADER3
for i = 1:11
    K(i)         = floor(1/nu(i));
    nuMod(i)     = (1-K(i)*nu(i))/1; 
    M(i,:)       = sqrt((1-nu(i)/2*(1+2*nu(i)-nu(i)^2)+nu(i).*cos(betaVector).*(2/3+nu(i)-2/3*nu(i)^2)+nu(i)/6.*cos(2.*betaVector).*(nu(i)^2-1)).^2+(1/6*nu(i).*sin(2.*betaVector).*(1-nu(i)^2)-nu(i)/3.*sin(betaVector).*(4-nu(i)^2)).^2);
    MMod(i,:)    = sqrt((1-nuMod(i)/2*(1+2*nuMod(i)-nuMod(i)^2)+nuMod(i).*cos(betaVector).*(2/3+nuMod(i)-2/3*nuMod(i)^2)+nuMod(i)/6.*cos(2.*betaVector).*(nuMod(i)^2-1)).^2+(1/6*nuMod(i).*sin(2.*betaVector).*(1-nuMod(i)^2)-nuMod(i)/3.*sin(betaVector).*(4-nuMod(i)^2)).^2);
    MTotal(i,:)  = M(i,:).^K(i).*MMod(i,:);
end
titulo = 'UWC3-ADER3';
tituloFigura = 'UWC3_DT-ADER3';
plotnotpolar(betaVector,MTotal,titulo,tituloFigura,'A',31);
plotnotpolar(betaVector,MTotal,titulo,tituloFigura,'Z',32);
% Representacion del número de pasos temporales necesarios
grafica = figure(33)
bar(nu,log10(K));
hold on;
xlabel('\boldmath$$\sigma$$','Interpreter','latex');
ylabel('\boldmath$$log_{10}(K)$$','Interpreter','latex');
xticks([0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0]);
xticklabels({'0.001' '0.1' '0.2' '0.3' '0.4' '0.5' '0.6' '0.7' '0.8' '0.9' '1.0'});
grid on;
for i1=1:numel(K)
    text(nu(i1),log10(K(i1)),num2str(K(i1),'%d'),...
               'HorizontalAlignment','center',...
               'VerticalAlignment','bottom')
end
saveas(grafica,strcat(tituloFigura, '_pasosTemporales'),'epsc')
hold off;
%% Estudio de cómo afecta el paso temporal - UCW3-RK3
for i = 1:11
    K(i)         = floor(1/nu(i));
    nuMod(i)     = (1-K(i)*nu(i))/1; 
    M(i,:)       = sqrt(((1-nu(i)/2-5*nu(i)^2/24+59*nu(i)^3/432) + (2*nu(i)/3-5*nu(i)^2/18-4*nu(i)^3/72).*cos(betaVector) + (-nu(i)/6+23*nu(i)^2/36-35*nu(i)^3/144).*cos(2.*betaVector) + (-nu(i)^2/6+155*nu(i)^3/648).*cos(3.*betaVector) + (nu(i)^2/72-13*nu(i)^3/144).*cos(4.*betaVector) + (nu(i)^3/72).*cos(5.*betaVector) - (nu(i)^3/1296).*cos(6.*betaVector)).^2+((-4*nu(i)/3+11*nu(i)^2/18+6*nu(i)^3/72).*sin(betaVector) + (nu(i)/6-19*nu(i)^2/36+27*nu(i)^3/144).*sin(2.*betaVector) + (nu(i)^2/6-163*nu(i)^3/648).*sin(3.*betaVector) + (-nu(i)^2/72+13*nu(i)^3/144).*sin(4.*betaVector) - (nu(i)^3/72).*sin(5.*betaVector) + (nu(i)^3/1296).*sin(6.*betaVector)).^2);
    MMod(i,:)    = sqrt(((1-nuMod(i)/2-5*nuMod(i)^2/24+59*nuMod(i)^3/432) + (2*nuMod(i)/3-5*nuMod(i)^2/18-4*nuMod(i)^3/72).*cos(betaVector) + (-nuMod(i)/6+23*nuMod(i)^2/36-35*nuMod(i)^3/144).*cos(2.*betaVector) + (-nuMod(i)^2/6+155*nuMod(i)^3/648).*cos(3.*betaVector) + (nuMod(i)^2/72-13*nuMod(i)^3/144).*cos(4.*betaVector) + (nuMod(i)^3/72).*cos(5.*betaVector) - (nuMod(i)^3/1296).*cos(6.*betaVector)).^2+((-4*nuMod(i)/3+11*nuMod(i)^2/18+6*nuMod(i)^3/72).*sin(betaVector) + (nuMod(i)/6-19*nuMod(i)^2/36+27*nuMod(i)^3/144).*sin(2.*betaVector) + (nuMod(i)^2/6-163*nuMod(i)^3/648).*sin(3.*betaVector) + (-nuMod(i)^2/72+13*nuMod(i)^3/144).*sin(4.*betaVector) - (nuMod(i)^3/72).*sin(5.*betaVector) + (nuMod(i)^3/1296).*sin(6.*betaVector)).^2);
    MTotal(i,:)  = M(i,:).^K(i).*MMod(i,:);
end
titulo = 'UWC3-RK3';
tituloFigura = 'UWC3_DT-RK3';
plotnotpolar(betaVector,MTotal,titulo,tituloFigura,'A',34);
plotnotpolar(betaVector,MTotal,titulo,tituloFigura,'Z',35);
% Representacion del número de pasos temporales necesarios
grafica = figure(36)
bar(nu,log10(K));
hold on;
xlabel('\boldmath$$\sigma$$','Interpreter','latex');
ylabel('\boldmath$$log_{10}(K)$$','Interpreter','latex');
xticks([0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0]);
xticklabels({'0.001' '0.1' '0.2' '0.3' '0.4' '0.5' '0.6' '0.7' '0.8' '0.9' '1.0'});
grid on;
for i1=1:numel(K)
    text(nu(i1),log10(K(i1)),num2str(K(i1),'%d'),...
               'HorizontalAlignment','center',...
               'VerticalAlignment','bottom')
end
saveas(grafica,strcat(tituloFigura, '_pasosTemporales'),'epsc')
hold off;