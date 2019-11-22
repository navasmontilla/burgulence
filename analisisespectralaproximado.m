clear all;close all;
CFL           = [0.001]; %[0.001 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0];
NCELLS        = 200;
NMAX          = NCELLS/2;
a             = 1.0;
% Lista de esquemas:
%   1:  FOU
%   2:  ADER3 - UWC3
%   3:  ADER5 - UWC5
%   4:  ADER7 - UWC7
%   5:  ADER9 - UWC9
%   6:  ADER3 - WENO3
%   7:  ADER5 - WENO5
%   8:  ADER7 - WENO7
%   9:  ADER9 - WENO9
%   10: RK3   - UWC3
%   11: RK3   - UWC5
%   12: RK3   - UWC7
%   13: RK3   - WENO3
%   14: RK3   - WENO5
%   15: RK3   - WENO7
%   16: FOU
scheme          = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16];
order           = [1 3 5 7 9 3 5 7 9  3  5  7  3  5  7  1];  % El orden de la discretización espacial
spatialScheme   = [1 3 5 7 9 3 5 7 9  1  1  1  2  2  2  0];  % El 1 y el 2 de las cuatro últimas posiciones se corresponden con los grados 3 y 5 del UWC y del WENO en el código semidsicreto.
reconstruction  = [4 4 4 4 4 1 1 1 1  0  0  0  0  0  0  0];  % Tipo de reconstruccion: 1-WENO SHU, 2-WENO PW, 3-WENO Z, 4-OPTIMAL WEIGHTS (UWC)
aderDegree      = [1 2 3 4 5 2 3 4 5  0  0  0  0  0  0  0];
TMAX            = 0.1;
L               = 1.0;
mallado         = 2;
NCELLY          = 10; % En principio no nos importa en 1d
color           = distinguishable_colors(length(CFL)); % Colores

for mm = 1:length(CFL)
   legendInfo1(mm,:) = sprintf('CFL=%.3f',CFL(mm));
   legendInfo2(mm,:) = sprintf('CFL=%.3f',CFL(mm));
end
legendInfo1(length(CFL)+1,:) = sprintf('--Ideal--');
legendInfo1 = cellstr(legendInfo1);
legendInfo2 = cellstr(legendInfo2);

for jj=10:10
    for mm = 1:length(CFL)
        if scheme(jj) < 9
            fichero1 = fopen('dataADER.txt','w');
            fprintf(fichero1,'0.1D0\t !TMAX\n');
            fprintf(fichero1,'%f\t !CFL\n',CFL(mm));
            fprintf(fichero1,'1.D0\t !LAMBDA X\n');
            fprintf(fichero1,'0.D0\t !LAMBDA Y\n');
            fclose(fichero1);
        end
        contadorDisp = 0;
        contadorDiff = 0;
        for n=1:NMAX
            % ESCRIBIR EL FICHERO
            fprintf('scheme = %d - CFL = %f - n = %f\n',jj,CFL(mm),n);
            digits(16);
            DX(n) = (2*pi*n/(2*NMAX));   % Paso temporal para cada n, 16 decimales
            DT(n) = DX(n)*CFL(mm)/a;
            if scheme(jj) < 9
                disp('ADER');
                fichero2 = fopen('dataWENO.txt','w');
                fprintf(fichero2,'%i\t !NCELLX\n',NCELLS);
                fprintf(fichero2,'%i\t !NCELLY\n',NCELLY);
                fprintf(fichero2,'%i\t !K >> el orden es 2*k-1\n',aderDegree(jj));
                fprintf(fichero2,'%.16f\t !DELTAX PONER 16 CIFRAS\n',DX(n));
                fprintf(fichero2,'1.D-25\t	!EPSILON nada\n');
                fprintf(fichero2,'%i\t !SELECTOR MALLADO: 1=CONSTANTE, 2=VARIABLE CON REFINAMIENTO EN EXTREMOS nada\n',mallado);
                fprintf(fichero2,'%i\t !TIPO DE reconstruccion: 1-WENO SHU, 2-WENO PW, 3-WENO Z, 4-OPTIMAL WEIGHTS (UWC)\n',reconstruction);
                fprintf(fichero2,'1.D-8\t !THRESHOLD FOR THETA COMPUTATION. DIFERENCIA ENTRE U(N+1) Y U(N) MINIMA CON LA QUE HACER THETA=0 nada\n');
                fclose(fichero2);
                system('ader.exe');
            else
                %disp('semidiscreto');
                fichero2 = fopen('data.txt','w');
                fprintf(fichero2,'%d\n',NCELLS);
                fprintf(fichero2,'%.16f\n',DX(n));
                fprintf(fichero2,'%f\n',a);
                fprintf(fichero2,'%f\n',CFL(mm));
                fprintf(fichero2,'%d\n',order(jj));
                fprintf(fichero2,'%d\n',spatialScheme(jj));
                fclose(fichero2);
                system('discreto7.exe');
            end

%             % IC writing
%             fichero3 = fopen('dataCI.txt','w');
%             fichero4 = fopen('dataEXACT.txt','w');
%             for ii = 1:NCELLS
%                 xcell(ii) = DX(n)*(ii) - DX(n)/2.0;
%                 x1        = xcell(ii)  + DX(n)/2.0;
%                 x2        = xcell(ii)  - DX(n)/2.0;
%                 ui(ii)    = (-cos(2*pi*x1/L) + cos(2*pi*x2/L))/(DX(n)*2*pi);
%                 fprintf(fichero3,'%.16f\n',ui(ii));
%                 % Exact solution
%                 digits(16);
%                 uf(ii) = (-cos(2*pi*(x1-a*TMAX)/L) + cos(2*pi*(x2-a*TMAX)/L))/(DX(n)*2*pi);
%                 fprintf(fichero4,'%.16f\n',uf(ii));
%             end
%             fclose(fichero3);
%             fclose(fichero4);

            % LECTURA DE RESULTADOS
            A  = dlmread('SOLUTION1dX.dat','%t');
            B  = dlmread('SOLUTION1dX_CI.dat','%t');
            x  = B(:,1);
            u  = A(:,2);
            ui = B(:,2);

            % Transformada rápida de Fourier de las CI
            L1  = length(ui);    % Length of signal
            ts1 = (0:L1-1)*DX(n);   % space vector
            Ui  = fft(ui);       % Fourier transform

            % Normalización de los datos
            Ui_bar         = zeros(L1/2+1,1);
            Ui_bar(2:L1/2) = 2*real(Ui(2:L1/2))/NCELLS + 1i*2*imag(Ui(2:L1/2))/NCELLS;
            Ui_bar(1)      = real(Ui(1))/NCELLS        + 1i*2*imag(Ui(1))/NCELLS;
            Ui_bar(L1/2+1) = real(Ui(L1/2+1))/NCELLS   + 1i*2*real(Ui(L1/2+1))/NCELLS;
            % Máxima componente
            idx  = find(abs(Ui_bar)>0.01);
            imin = min(idx);

            % Transformada rápida de Fourier
            U = fft(u);        % Fourier transform
            % Normalización de los datos
            U_bar         = zeros(L1/2+1,1);
            U_bar(2:L1/2) = 2*real(U(2:L1/2))/NCELLS + 1i*2*imag(U(2:L1/2))/NCELLS;
            U_bar(1)      = real(U(1))/NCELLS        + 1i*2*imag(U(1))/NCELLS;
            U_bar(L1/2+1) = real(U(L1/2+1))/NCELLS   + 1i*2*real(U(L1/2+1))/NCELLS;

            % Cálculo de la longitud de onda modificada
            phi(n)   = 1i*log(U_bar(imin)/Ui_bar(imin))/(CFL(mm));
            xi(n)    = real(phi(n));
            eta(n)   = imag(phi(n));
            M_exp(n) = exp(-1i*phi(n)*CFL(mm));
            mag(n)   = abs(M_exp(n));
             % Comprobamos que se obtiene la misma frecuencia
            frecuencias = sprintf('kDx = %f - DFT = %f',DX(n),2*real(Ui_bar(imin)));
            disp(frecuencias);
            beta(n)= 2*pi*n/(2*NMAX);
            
            if abs(xi(n) - beta(n))/beta(n) > 0.01 & contadorDisp == 0
                porcientoDisp(n,jj,mm) = xi(n);
                contadorDisp = contadorDisp + 1;
            end
			if abs(1.0 - mag(n))/1.0 > 0.01 & contadorDiff == 0
                porcientoDiff(n,jj,mm) = xi(n);
                contadorDiff = contadorDiff + 1;
            end

        end
        % Representación de dispersión y difusión semidiscreta

        ct1 	= 0;
        xi_fix 	= xi;
        mag_fix	= mag;
        eta_fix = eta;
        der 	= zeros(1,length(xi));
        der2 	= zeros(1,length(xi));
        der3    = zeros(1,length(xi));
        tol 	= 0.3; %tolerancia
        for kk=2:1:length(xi)-1
            der(kk) 	= xi(kk) - xi(kk-1);
            der2(kk-1) 	= der(kk) - der(kk-1);
        end
        for kk=2:1:length(xi)-1
            dif_dos = abs(der2(kk-1)-der2(kk+1)); %diferencia entre picos a los lados
            sgn_dos = der2(kk-1)*der2(kk+1); %es >0 si los picos tienen mismo signo
            sgn1 	= der2(kk)*der2(kk+1); %sera -1 si es el pico central
            sgn2 	= der2(kk)*der2(kk-1); %sera -1 si es el pico central
            der2(kk-1);
            if dif_dos<tol*abs(der2(kk-1)) && sgn_dos>0 && sgn1<0 && sgn2<0 %que la diferencia entre picos a los lados sea pequeña
                xi_fix (kk) = 0.5*(xi(kk+1)  +  xi(kk-1) );
                mag_fix(kk) = 0.5*(mag(kk+1) +  mag(kk-1));
                eta_fix(kk) = 0.5*(eta(kk+1) +  eta(kk-1));
                ct1 	 	= ct1 + 1;
            end
        end

        % filter kernel
        G_num  = phi./beta;
        G_num2 = (xi+(eta)*i)./beta;

        ct1=0;
        im=imag(phi);
        im_fix=im;
        der=zeros(1,length(im));
        der2=zeros(1,length(im));
        tol=0.2; %tolerancia
        for kk=2:1:length(im)-1
            der(kk)=im(kk)-im(kk-1);
            der2(kk-1)=der(kk)-der(kk-1);
        end
        for kk=2:1:length(im)-1
            dif_dos=abs(der2(kk-1)-der2(kk+1)); %diferencia entre picos a los lados
            sgn_dos=der2(kk-1)*der2(kk+1); %es >0 si los picos tienen mismo signo
            sgn1=der2(kk)*der2(kk+1); %sera -1 si es el pico central
            sgn2=der2(kk)*der2(kk-1); %sera -1 si es el pico central
            der2(kk-1);
            if dif_dos<tol*abs(der2(kk-1)) && sgn_dos>0 && sgn1<0 && sgn2<0 %que la diferencia entre picos a los lados sea pequeña
                im_fix (kk) =0.5*(im(kk+1)  +  im(kk-1) );
                ct1=ct1+1;
            end
        end

        G_num_fix=(xi_fix+(im_fix)*1i)./beta; %ojo si has añadido el "cero", aquí no se tiene en cuenta.
        G_exact_1st=(sin(beta)+1i*(cos(beta)-1))./beta;

        nc=100;
        beta_long     =2*pi/(2*NMAX)*(0:1:((2*NMAX)*nc+NMAX)); %vector muy largo de betas para luego integrar
        full_curve    =zeros(1,length(beta_long));             %curva de phi en ese vector
        G_num_fix_long=zeros(1,length(beta_long));

        beta_single =zeros(1,length(beta)*2+1);             %vector simetrico de betas alrededor de cero [-pi,pi]. Aqui se incluye el cero
        single_curve=zeros(1,length(beta)*2+1);             %curva de phi en [-pi,pi]
        single_curve(NMAX+1)=0+0*1i;
        beta_single(NMAX+1)=0.0;
        for kk=1:length(beta)
            single_curve(NMAX+1+kk)=xi_fix(kk)+(im_fix(kk))*1i;
            single_curve(NMAX+1-kk)=xi_fix(kk)+(im_fix(kk))*1i;
            beta_single(NMAX+1+kk)=beta(kk);
            beta_single(NMAX+1-kk)=-beta(kk);
        end

        full_curve(1:NMAX+1)=single_curve(NMAX+1:end); %aqui se rellena la curva larga considerando periodicidad par en n*pi
        for r=1:nc
            ini=NMAX+1+(r-1)*(2*NMAX);
            fin=ini+2*NMAX;
            full_curve(ini:fin)=single_curve(1:end);
        end

        G_long_abs_sq=(abs(full_curve./beta_long)).^2;
        G_long_abs_sq(1)=1;

        %calculo de la integral numerica con regla trapecio de un semieje
        G_int=0.0;
        for kk=2:1:length(beta_long)
            aux1=0.5*(G_long_abs_sq(kk)+G_long_abs_sq(kk-1));
            aux2=(beta_long(kk)-beta_long(kk-1));
            G_int=G_int+aux1*aux2;
        end
        G_int2=2*G_int;
        if mm == 1
            rCoef(jj)=(2*pi/G_int2); %este es el indicador de la pagina 78 del paper (2pi en la fórmula va dividiendo!!)
        end
            
        % Añadir el punto (0,0)
        beta_plot(1)                      = 0.0;
        beta_plot(2:length(beta)+1)       = beta(1:length(beta));
        xi_plot(1)                        = 0.0;
        xi_plot(2:length(xi)+1)           = xi(1:length(xi));
        xi_fix_plot(1)                    = 0.0;
        xi_fix_plot(2:length(xi_fix)+1)   = xi_fix(1:length(xi_fix));
        mag_fix_plot(1)                   = 1.0;
        mag_fix_plot(2:length(mag_fix)+1) = mag_fix(1:length(mag_fix));
        eta_plot(1)                       = 0.0;
        eta_plot(2:length(eta)+1)         = eta(1:length(eta));
        eta_fix_plot(1)                   = 0.0;
        eta_fix_plot(2:length(eta_fix)+1) = eta_fix(1:length(eta_fix))./beta_plot(2:length(beta_plot));
        G_num_plot(1)                     = 1.0;
        G_num_plot(2:length(G_num_fix)+1) = G_num_fix(1:length(G_num_fix));

        % Representaciones
        figure3 = figure(3);
        p1=plot(beta_plot,xi_fix_plot,'Color',color(mm,:),'LineWidth',1.5);
        % set(get(get(p1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
        hold on;
        xlabel('\boldmath$$k\Delta x$$','Interpreter','latex');
        ylabel('\boldmath$$\tilde{k}\Delta x$$','Interpreter','latex');
        grid on;
        xlim([0 pi]);
        ylim([0 2.5]);
        set(gca,'XTick',0:pi/4:pi);
        set(gca,'XTickLabel',{'0','\pi/4','\pi/2','3\pi/4','\pi'});
        

        figure4 = figure(4);
        p2 = plot(beta_plot,abs(mag_fix_plot),'LineWidth',1.5,'Color',color(mm,:));
        % set(get(get(p2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
        hold on;
        xlim([0 pi]);
        ylim([0 1.2]);
        grid on;
        ylabel('\boldmath$$|M|$$','Interpreter','latex');
        xlabel('\boldmath$$k\Delta x$$','Interpreter','latex');
        set(gca,'XTick',0:pi/4:pi);
        set(gca,'XTickLabel',{'0','\pi/4','\pi/2','3\pi/4','\pi'});

        figure5 = figure(5);
        p3 = plot(beta_plot,abs(G_num_plot).^2,'LineWidth',1.5,'Color',color(mm,:));
        % set(get(get(p3,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
        hold on;
        grid on;
        xlabel('\boldmath$$k\Delta x$$','Interpreter','latex');
        ylabel('\boldmath$$|G|$$','Interpreter','latex');
        xlim([0 pi]);
        ylim([0.0 1.2]);
        set(gca,'XTick',0:pi/4:pi);
        set(gca,'XTickLabel',{'0','\pi/4','\pi/2','3\pi/4','\pi'});

    end
    
    porcientoDisp(porcientoDisp == 0) = nan;
    porcientoDiff(porcientoDiff == 0) = nan;
    for mm = 1:length(CFL)
        [longOndaDisp(jj,mm),indexDisp(jj,mm)] = min(porcientoDisp(:,jj,mm));
        diffOver1(jj,mm)                       = eta_fix_plot(indexDisp(jj,mm));
        errorDiff(jj,mm) = abs(1.0 - mag(indexDisp(jj,mm)))/1.0*100;
        [longOndaDiff(jj,mm),indexDiff(jj,mm)] = min(porcientoDiff(:,jj,mm));
        dispOver1(jj,mm)                       = xi_fix_plot(indexDiff(jj,mm));
        errorDisp(jj,mm) = abs(dispOver1(jj,mm) - beta(indexDiff(jj,mm)))/beta(indexDiff(jj,mm))*100;
    end
    

    figure(3)
    line([0 2*pi],[0 2*pi],'LineStyle',':','Color','black','LineWidth',0.5);
    lgd = legend('0.001','0.100','0.200','0.300','0.400','0.500','0.600','0.700','0.800','0.900','1.000','Ideal','Location','southwest');
    lgd.Title.String = 'CFL';
    hold off;
    
    figure(4)
    lgd = legend('0.001','0.100','0.200','0.300','0.400','0.500','0.600','0.700','0.800','0.900','1.000','Location','southwest');
    lgd.Title.String = 'CFL';
    hold off;
    
    figure(5);
    lgd = legend('0.001','0.100','0.200','0.300','0.400','0.500','0.600','0.700','0.800','0.900','1.000','Location','southwest');
    lgd.Title.String = 'CFL';
    hold off;

    saveas(figure3,sprintf('approx_dispersion_scheme_%i',jj),'epsc');
    saveas(figure4,sprintf('approx_diffusion_scheme_%i',jj),'epsc');
    saveas(figure5,sprintf('approx_G_num_scheme_%i',jj),'epsc');

end
