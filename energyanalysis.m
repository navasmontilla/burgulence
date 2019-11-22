clc; clear all; close all;

NFILES = 149;
START  = 50;
DT     = 4;
pend1  = -5.0/3.0;
pend2  = -2.0;
N1porc = 0.0;
DK     = 0.0;
k1disp = 0.0;
k1diff = 0.0;
NCELLS  = 4096;
Nc      = 80;
k1disp  = 0.879645943005142;
k1diff  = 0.879645943005142;
x_line1 = 1.0;
y_line1 = -3.0;
x_line2 = 2.75;
x_pend1 = 1.5;
y_pend1 = -3.2;
x_pend2 = 2.25;
y_pend2 = -4.5;
xlim_l  = 0.5;
xlim_r  = 3.5;
ylim_l  = -14.0;
ylim_r  = 0.0;
xlimf_l = 0.5;
xlimf_r = 3.5;
ylimf_l = -14.0;
ylimf_r = 0.0;
xzoom_l = log10(Nc);
xzoom_r = 3.3;
yzoom_b = -4.0;
yzoom_t = -1.0;
DX      = 1/NCELLS;
N1disp  = k1disp*NCELLS/(pi*2);
N1diff  = k1diff*NCELLS/(pi*2);

% Data reading
A    = dlmread(sprintf('SOLUTION1dX_CI.dat'),'%t');
x    = A(:,1);
u_CI = A(:,2);
for ii = 1:NFILES
     filename = sprintf('SOLUTION1dX_%d.txt',ii-1);
     B        = dlmread(filename,'%t');
     u(:,ii)  = B(:,2);
end

% Data normalization
U      = abs(fft(u(:,START:NFILES)));
U_norm = zeros(NCELLS/2+1,(NFILES-START+1));
for ii = 1:NFILES-START+1
    U_norm(2:NCELLS/2,ii) = 2*real(U(2:NCELLS/2,ii))/NCELLS - 1i*2*imag(U(2:NCELLS/2,ii))/NCELLS;
    U_norm(1,ii)          = real(U(1,ii))/NCELLS            - 1i*2*imag(U(1,ii))/NCELLS;
    U_norm(NCELLS/2+1,ii) = real(U(NCELLS/2+1,ii))/NCELLS   - 1i*2*imag(U(NCELLS/2+1,ii))/NCELLS;
end
E         = U_norm.^2/2.0;
Epromedio = sum(E,2)/(NFILES-START+1);
ts1       = (0:length(Epromedio)-1)*DX; % space vector
x_plot    = 1:length(Epromedio);

% Finding the slope
x_log       = log10(x_plot)';
E_log       = log10(abs(Epromedio));
leftlimit1  = x_line1;
rightlimit1 = log10(Nc);
leftlimit2  = log10(Nc);
rightlimit2 = x_line2; % Estos limites estan medido sobre el articulo Moura
index1      = find((x_log > leftlimit1) & (x_log < rightlimit1));
index2      = find((x_log > leftlimit2) & (x_log < rightlimit2));
slope1      = polyfit(x_log(index1),E_log(index1),1);
E_poly1     = polyval(slope1,x_log(index1));
slope2      = polyfit(x_log(index2),E_log(index2),1);
E_poly2     = polyval(slope2,x_log(index2));

% Finding the relative error

[minValueDisp,closestIndexDisp] = min(abs(x_log - log10(N1disp)));
closestValueDisp = x_log(closestIndexDisp);
[minValueDiff,closestIndexDiff] = min(abs(x_log - log10(N1diff)));
closestValueDiff = x_log(closestIndexDiff);
E_disp = polyval(slope2,log10(N1disp));
E_diff = polyval(slope2,log10(N1diff));
y_disp = E_log(closestIndexDisp);
y_diff = E_log(closestIndexDiff);
errorDisp = abs((y_disp - E_disp)/E_disp*100);
errorDiff = abs((y_diff - E_diff)/E_diff*100);

% Setting turbulent state
tfluctuation = 0:DT:(DT*NFILES);
fluctuation1 = abs(u_CI -1.0);
fluctuation  = abs(u - 1.0);
fluctProm(1) = sum(fluctuation1)/NCELLS;
fluctProm(2:NFILES+1) = sum(fluctuation)/NCELLS;

% Results writing
fichero2 = fopen(sprintf('energycascade.txt'),'w');
fichero3 = fopen(sprintf('energycascadelog.txt'),'w');
fichero4 = fopen(sprintf('energy.txt'),'w');
for ii = 1:length(x_plot)
    fprintf(fichero2,'%.16f %.16f\n',x_plot(1,ii), abs(Epromedio(ii)));
    fprintf(fichero3,'%.16f %.16f\n',x_log(ii),E_log(ii));
    fprintf(fichero4,'%.16f\n',Epromedio(ii));
end
fclose(fichero2);
fclose(fichero3);
fclose(fichero4);

% Plotting

figure1 = figure(1);
plot(tfluctuation,fluctProm,'LineWidth',1.5);
hold on;
grid on;
ylim([0.0 0.16]);
set(gca,'YTick',0:0.04:0.16);
xlabel('\boldmath$$t$$','Interpreter','latex');
ylabel('\boldmath$$u(0,t)^,$$','Interpreter','latex');
hold off;
saveas(figure1,sprintf('fluctuations'),'png');
saveas(figure1,sprintf('fluctuations'),'epsc');

figure2 = figure(2);
plot(x,u_CI,'Color',[0 0 1],'LineWidth',1.2);
hold on;
%ylim([0.8 1.2]);
plot(x,u(:,NFILES-1),'Color','black','LineWidth',1.2);
hold off;
ylim([0.7 1.3]);
set(gca,'YTick',0.7:0.15:1.3);
set(gcf, 'units', 'normalized');
set(gcf, 'Position', [0.1, 0.1, 0.4, 0.2]);
saveas(figure2,sprintf('final'),'png');
saveas(figure2,sprintf('final'),'epsc');

filtrado = sgolayfilt(u(:,NFILES-1),1,1201);  


figure11 = figure(11);
plot(x,u_CI,'Color',[0 0 1],'LineWidth',1.2);
hold on;
%ylim([0.8 1.2]);
plot(x,u(:,NFILES-1),'Color','black','LineWidth',1.2);
plot(x,filtrado,'Color','red','LineWidth',1.2);
hold off;
ylim([0.75 1.25]);
set(gca,'YTick',0.75:0.125:1.25);
xlabel('\boldmath$$x$$','Interpreter','latex');
ylabel('\boldmath$$u(x)$$','Interpreter','latex');
set(gcf, 'units', 'normalized');
set(gcf, 'Position', [0.1, 0.1, 0.4, 0.2]);
saveas(figure11,sprintf('finalfiltrado'),'png');
saveas(figure11,sprintf('finalfiltrado'),'epsc');

figure10 = figure(10);
plot(x,u_CI,'Color',[0 0 1],'LineWidth',1.2);
hold on;
%ylim([0.8 1.2]);
plot(x,u(:,NFILES-1),'Color','black','LineWidth',1.2);
hold off;
xlim([0.0 0.4]);
ylim([0.9 1.15]);
set(gca,'YTick',0.7:0.15:1.3);
set(gcf, 'units', 'normalized');
set(gcf, 'Position', [0.1, 0.1, 0.4, 0.2]);
saveas(figure10,sprintf('finalzoom'),'png');
saveas(figure10,sprintf('finalzoom'),'epsc');

upromedio = sum(u(:,START:NFILES),2)/(NFILES - START + 1);
figure3 = figure(3);
plot(x,u_CI);
hold on;
plot(x,upromedio,'--','Color','red');
hold off;
saveas(figure3,sprintf('promedio'),'png');
saveas(figure3,sprintf('promedio'),'epsc');

figure4 = figure(4);
plot(x_log,E_log,'Color','black');
hold on;
grid on;
%l1 = line([2.0 3.5],[-2 -2-5/3.0*(3.5 - 2.0)],'Color','black','LineWidth',1.5);
text(1.65,-2.4,'\boldmath$$-1.91$$','Interpreter','latex','Rotation',-15,'FontSize',12);
%l2 = line([1.5 3.0],[-6 -6-2.0*(3.0 - 1.5)],'Color','black','LineWidth',1.5);
text(1.8,-7.4,'\boldmath$$-2.10$$','Interpreter','latex','Rotation',-17,'FontSize',12);
%xlim([1.5 3.7]);
%xticks(1.5:0.5:3.5);
%ylim([-11 -4]);
%yticks(-10:5:-5);
plot(x_log(index1),E_poly1,'LineWidth',1.4);
plot(x_log(index2),E_poly2,'LineWidth',1.4);
plot(log10(N1disp),E_disp,'r*');
plot(log10(N1diff),E_diff,'r*','Color','blue');
line([log10(Nc) log10(Nc)],[-10 4],'LineStyle','--','Color','black');
hold off;
saveas(figure4,sprintf('cascada1'),'png');
saveas(figure4,sprintf('cascada1'),'epsc');

figure5 = figure(5);
plot(x_log,E_log,'Color','black');
hold on;
plot(x_log,E_log + (5/3).*x_log,'Color','blue');
plot(x_log,E_log + (2).*x_log,'Color','red');
grid on;
line([log10(Nc) log10(Nc)],[-10 4],'LineStyle','--','Color','black');
hold off;
saveas(figure5,sprintf('cascadaPendiente'),'png');
saveas(figure5,sprintf('cascadaPendiente'),'epsc');


figure6 = figure(6);
p61 = plot(x_log(1:length(x_log)-2),E_log(1:length(x_log)-2),'Color','black','LineWidth',1.1);
set(get(get(p61,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
hold on;
grid on;
xlim([xlim_l xlim_r]);
ylim([ylim_l ylim_r]);
set(gca,'YTick',ylim_l:2:0);
xlabel('\boldmath$$\log_{10}k$$','Interpreter','latex');
ylabel('\boldmath$$\log_{10}E(k)$$','Interpreter','latex');
line6_Nc   = line([log10(Nc) log10(Nc)],[ylim_l ylim_r],'LineStyle','--','Color','black');
line6_disp = line([log10(N1disp) log10(N1disp)],[ylim_l ylim_r],'LineStyle','-.','Color','black');
line6_diff = line([log10(N1diff) log10(N1diff)],[ylim_l ylim_r],'LineStyle',':','Color','black');
yyaxis right;
p62 = plot(x_log(1:length(x_log)-2),E_log(1:length(x_log)-2) + (2).*x_log(1:length(x_log)-2),'Color','red');
set(get(get(p62,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
ylim([ylim_l ylim_r]);
set(gca,'YColor',[1 0 0],'YTick',ylim_l:2:0);
ylabel('\boldmath$$\log_{10}(E(k)/k^{-2})$$','Interpreter','latex');
l61 = line([x_line1 log10(Nc)],[y_line1 y_line1+pend1*(log10(Nc) - x_line1)],'Color',[0.7 0.7 0.7],'LineWidth',1.5);
set(get(get(l61,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
text(x_pend1,y_pend1,'\boldmath$$k^{-5/3}$$','Interpreter','latex','Rotation',-13,'FontSize',12);
l62 = line([log10(Nc) x_line2],[y_line1+pend1*(log10(Nc) - x_line1) y_line1+pend1*(log10(Nc) - x_line1)+pend2*(x_line2 - log10(Nc))],'Color',[0.4 0.4 0.4],'LineWidth',1.5);
set(get(get(l62,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
text(x_pend2,y_pend2,'\boldmath$$k^{-2}$$','Interpreter','latex','Rotation',-15,'FontSize',12);
leg6 = legend('\boldmath$N_c$','\boldmath$k_{1\%}^{disp}$','\boldmath$k_{1\%}^{diff}$');
set(leg6,'Interpreter','latex','Location','southwest','FontSize',10);
saveas(figure6,sprintf('cascade'),'png');
saveas(figure6,sprintf('cascade'),'epsc');

figure7 = figure(7);
p71 = plot(x_log(1:length(x_log)-2),E_log(1:length(x_log)-2),'Color','black','LineWidth',1.3);
set(get(get(p71,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
hold on;
grid on;
xlim([xlimf_l xlimf_r]);
ylim([ylimf_l ylimf_r]);
set(gca,'YTick',ylim_l:2:0);
xlabel('\boldmath$$\log_{10}k$$','Interpreter','latex');
ylabel('\boldmath$$\log_{10}E(k)$$','Interpreter','latex');
line7_Nc   = line([log10(Nc) log10(Nc)],[ylimf_l ylimf_r],'LineStyle','--','Color','black','LineWidth',1.1);
line7_disp = line([log10(N1disp) log10(N1disp)],[ylimf_l ylimf_r],'LineStyle','-.','Color','black','LineWidth',1.1);
line7_diff = line([log10(N1diff) log10(N1diff)],[ylimf_l ylimf_r],'LineStyle',':','Color','black','LineWidth',1.1);
yyaxis right;
p72 = plot(x_log(1:length(x_log)-2),E_log(1:length(x_log)-2) + (2).*x_log(1:length(x_log)-2),'Color','red','LineWidth',1.3);
set(get(get(p72,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
ylim([ylimf_l ylimf_r]);
set(gca,'YColor',[1 0 0],'YTick',ylim_l:2:0);
ylabel('\boldmath$$\log_{10}(E(k)/k^{-2})$$','Interpreter','latex');
l71 = line([x_line1 log10(Nc)],[y_line1 y_line1+pend1*(log10(Nc) - x_line1)],'Color',[0.7 0.7 0.7],'LineWidth',1.3);
set(get(get(l71,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
text(x_pend1,y_pend1,'\boldmath$$k^{-5/3}$$','Interpreter','latex','Rotation',-13,'FontSize',12);
l72 = line([log10(Nc) x_line2],[y_line1+pend1*(log10(Nc) - x_line1) y_line1+pend1*(log10(Nc) - x_line1)+pend2*(x_line2 - log10(Nc))],'Color',[0.4 0.4 0.4],'LineWidth',1.3);
set(get(get(l72,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
text(x_pend2,y_pend2,'\boldmath$$k^{-2}$$','Interpreter','latex','Rotation',-15,'FontSize',12);
leg7 = legend('\boldmath$N_c$','\boldmath$k_{1\%}^{disp}$','\boldmath$k_{1\%}^{diff}$');
set(leg7,'Interpreter','latex','Location','southwest','FontSize',10);
set(gcf, 'units', 'normalized');
set(gcf, 'Position', [0.1, 0.1, 0.4, 0.2]);
saveas(figure7,sprintf('cascadeFlattened'),'png');
saveas(figure7,sprintf('cascadeFlattened),'epsc');

figure8 = figure(8);
plot(x_log,E_log + (2).*x_log,'Color','red');
hold on;
%l1 = line([log10(Nc) 3.5],[y0 y0],'Color',[0.7 0.7 0.7],'LineWidth',1.5);
grid on;
xlim([xzoom_l xzoom_r]);
ylim([yzoom_b yzoom_t]);
line([log10(N1porc) log10(N1porc)],[yzoom_b yzoom_t],'LineStyle',':','Color','black');
%set(gca,'YTick',1.0:0.25:2.25);
xlabel('\boldmath$$\log_{10}k$$','Interpreter','latex');
ylabel('\boldmath$$\log_{10}(E(k)/k^{-2})$$','Interpreter','latex');
saveas(figure8,sprintf('cascadezoom_%s'),'png');
saveas(figure8,sprintf('cascadezoom_%s'),'epsc');

figure9 = figure(9);
plot(x_log,E_log + (2).*x_log,'Color','red');
hold on;
%l1 = line([log10(Nc) 3.5],[y0+(2)*log10(Nc) y0+(2)*log10(Nc)],'Color',[0.7 0.7 0.7],'LineWidth',1.5);
grid on;
xlim([xzoom_l xzoom_r]);
ylim([yzoom_b yzoom_t]);
line([log10(N1porc) log10(N1porc)],[yzoom_b yzoom_t],'LineStyle',':','Color','black');
%set(gca,'YTick',1.0:0.25:2.25);
xlabel('\boldmath$$\log_{10}k$$','Interpreter','latex');
ylabel('\boldmath$$\log_{10}(E(k)/k^{-2})$$','Interpreter','latex');
set(gcf, 'units', 'normalized');
set(gcf, 'Position', [0.1, 0.1, 0.6, 0.3]);
saveas(figure9,sprintf('cascadezoomFlattened'),'png');
saveas(figure9,sprintf('cascadezoomFlattened'),'epsc');
