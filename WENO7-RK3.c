#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#define NTHREADS 23
#define tolerancia 0.0000000000001
#define openmp

/*******************************************************/
/***************** Function Definitions ****************/
/*******************************************************/

void rightLeft3(double *vector, int R_L, int R_C, int R_R, int L_L, int L_C, int L_R, double *out_R, double *out_L, int N, int spatialScheme);
void rightLeft5(double *vector, int R_LL, int R_L, int R_C, int R_R, int R_RR, int L_LL, int L_L, int L_C, int L_R, int L_RR, double *out_R, double *out_L, int N, int spatialScheme);
void rightLeft7(double *vector, int R_LLL, int R_LL, int R_L, int R_C, int R_R, int R_RR, int R_RRR, int L_LLL, int L_LL, int L_L, int L_C, int L_R, int L_RR, int L_RRR, double *out_R, double *out_L, int N, int spatialScheme);
void fluxCalculus(double *vector, double *in_R, double *in_L, double *out_R, double *out_L, int N, double *lambda, int equation);
void weno5(double *vector, int R_LL, int R_L, int R_C, int R_R, int R_RR, int L_LL, int L_L, int L_C, int L_R, int L_RR, double *out_R, double *out_L, int N, int spatialScheme);

/*******************************************************/
/********************* Main Program ********************/
/*******************************************************/

int main(){
// Definition of read variables
double TMAX, CFL, DX, t, DT, a, a_plus, a_minus, control[3], burgersControl[3], L_dom;
int NCELL, i, k, degree, equation, sourceTerm, schemeOrder, spatialOrder, spatialScheme, Nc;
FILE *CI, *SF, *SF_val, *data, *dataCI;

// Inicialization
TMAX          = 0.0;
CFL           = 0.0;
NCELL         = 0;
DX            = 0.0;
t             = 0.0;
DT            = 0.0;
a             = 0.0;
a_plus        = 0.0;
a_minus       = 0.0;
degree        = 0;
equation      = 0;
sourceTerm    = 0;
schemeOrder   = 0;
spatialOrder  = 0;
spatialScheme = 0;
Nc 			  = 0;
burgersControl[0] = 0.0;
burgersControl[1] = 1.0;
for(i=0;i<3;i++){
    control[i]        = 0.0;
    burgersControl[i] = 0.0;
}

#ifdef openmp
	omp_set_num_threads(NTHREADS);
	printf("The number of threads is set to %d.\n",NTHREADS);
#endif

// Data reading
data = fopen("data.txt", "r");
if (data == NULL){
      perror("Error al abrir el archivo data.txt.\n");
      exit(EXIT_FAILURE);
}
fscanf(data,"%d",&NCELL);
fscanf(data,"%lf",&CFL);
fscanf(data,"%lf",&TMAX);
fscanf(data,"%d",&equation);      // Type of equation              -> 1: advective, 2: Burgers
fscanf(data,"%d",&sourceTerm);    // Source term switch            -> 0: No, 1: Yes
fscanf(data,"%d",&schemeOrder);   // Time discretization order     -> 1: FOU, 2: RK2, 3: RK3
fscanf(data,"%d",&spatialOrder);  // Spatial discretization order  -> 3: UWC3/WEN03, 5: UWC5/WENO5
fscanf(data,"%d",&spatialScheme); // Spatial discretization scheme -> 1: UWC, 2: WENO
fscanf(data,"%d",&Nc); // Spatial discretization scheme -> 1: UWC, 2: WENO
fclose(data);
printf("NCELL = %d\n",NCELL);
printf("source = %d\n",sourceTerm);
printf("Nc = %d\n",Nc);
// Data values
L_dom = 2.0;
DX = L_dom/NCELL;
printf("DX = %f\n",DX);
a_plus  = (a+fabs(a))/2.0;
a_minus = (a-fabs(a))/2.0;
degree  = 3;
burgersControl[2] = 1.0;

/************* Initial conditions ************/
/*
dataCI = fopen("dataCI.txt","r");
if (dataCI == NULL){
      perror("Error al abrir el archivo data.txt.\n");
      exit(EXIT_FAILURE);
}
double u[NCELL];
for(i=0;i<NCELL;i++)
    u[i] = 0.0;
for(i=0;i<NCELL;i++){
    fscanf(dataCI,"%lf",&u[i]);
}
fclose(dataCI);
// IC writing
CI = fopen("SOLUTION1dX_CI.dat", "w");
for(i=0;i<NCELL;i++){
    fprintf(CI,"%d\t %.16f\n",i,u[i]);
}
fclose(CI);
*/

double xcell[NCELL], u[NCELL], x1, x2;
x1            = 0.0;
x2            = 0.0;
for(i=0;i<NCELL;i++)
    xcell[i] = 0.0;
/*
CI = fopen("SOLUTION1dX_CI.dat", "w");
for(i=0;i<NCELL;i++){
    xcell[i] = DX/2*(2*i+1);
    if(i < 7 || (i > 25 && i < 50) ){
        u[i] = 2.0;
    }else{
        u[i] = 0.0;
    }
    fprintf(CI,"%.16f\t %.16f\n",xcell[i],u[i]);
}
fclose(CI);
*/
/*
CI = fopen("SOLUTION1dX_CI.dat", "w");
for(i=0;i<NCELL;i++){
    xcell[i] = DX/2*(2*i+1);
    x1   = xcell[i] + DX/2.0;
    x2   = xcell[i] - DX/2.0;
    u[i] = 1 + (-cos(2*M_PI*x1)+cos(2*M_PI*x2))/(DX*2*M_PI);
    fprintf(CI,"%.16f\t %.16f\n",xcell[i],u[i]);
}
fclose(CI);
*/

CI = fopen("SOLUTION1dX_CI.dat", "w");
for(i=0;i<NCELL;i++){
    xcell[i] = -1 + DX/2*(2*i+1);
    u[i] = 1.0;
    fprintf(CI,"%.16f\t %.16f\n",xcell[i],u[i]);
	//printf("%.16f\t %.16f\n",xcell[i],u[i]);
}
fclose(CI);

/****************** CFL Condition ****************/
double uDT;
if(equation == 1){
    DT = DX*CFL/fabs(a);
}

// More definitions
double unew[NCELL], uval[NCELL], unew_val[NCELL], u_R[NCELL], u_L[NCELL];
double f_L[NCELL], f_R[NCELL], du[NCELL], flux[NCELL], u_mean[NCELL], u_plus[NCELL], u_minus[NCELL];
double lambda[NCELL], u1[NCELL], u2[NCELL], newDT;
int step[schemeOrder];
newDT = 0.0;
for(i=0;i<NCELL;i++){
    unew[i]  	= 0.0;
	uval[i] 	= 0.0;
	unew_val[i] = 0.0;
	u_R[i]      = 0.0;
	u_L[i]      = 0.0;
	f_R[i]      = 0.0;
	f_L[i]      = 0.0;
	flux[i]     = 0.0;
	u_mean[i]   = 0.0;
	u_plus[i]   = 0.0;
	u_minus[i]  = 0.0;
	lambda[i]   = 0.0;
	u1[i]       = 0.0;
	u2[i]       = 0.0;
}
for(i=0;i<NCELL;i++){
	uval[i]   = u[i];
	lambda[i] = a;
}
for(i=0;i<schemeOrder;i++){
    step[i] = i + 1;
}
// Actualization equation coefficients
long double a_eq[3], b_eq[3], c_eq[3];
for(i=0;i<3;i++){
    a_eq[i] = 0.0;
    b_eq[i] = 0.0;
    c_eq[i] = 0.0;
}
a_eq[0] = 1.0;
a_eq[1] = 3.0/4.0;
a_eq[2] = 1.0/3.0;
b_eq[0] = 0.0;
b_eq[1] = 1.0/4.0;
b_eq[2] = 1.0 - a_eq[2]; // Necessary defining fractions this way
c_eq[0] = 1.0;
c_eq[1] = 1.0/4.0;
c_eq[2] = 1.0 - a_eq[2];

// Source term

int n_source, jj;
double amp_source, a_source, b_source, source[NCELL], suma_source, L_source, den_aux, x1_source, x2_source;
amp_source  = 0.04;
L_source    = 2.0;
a_source    = 0.0;
b_source    = 0.0;
x1_source   = 0.0;
x2_source   = 0.0;
suma_source = 0.0;
den_aux     = 0.0;
for(i=0;i<NCELL;i++){
    source[i] = 0.0;
}
double tdumping, taux, tthreshold;
int nt, ntwrite;
nt 		   = 0;
ntwrite    = 0;
taux 	   = 0.0;
tdumping   = 0.0;
tthreshold = 0.0;
tdumping   = 4.0;
tthreshold = 0.0;
char filename[64];

// Random array reading
double *randomarray, aux_array;
int count_array, point_array, N_random, index_rand;
FILE *randArray;
aux_array   = 0.0;
count_array = 0;
point_array = 0;
N_random 	= 6000000;
index_rand  = 0;
printf("Antes del fichero\n");
randomarray = (double*)malloc(N_random*sizeof(double));
randArray   = fopen("randomarray.txt","r");
while (fscanf(randArray, "%lf", &aux_array)!= EOF){
	randomarray[count_array] = aux_array;
	count_array    = count_array + 1;
}
fclose(randArray);
printf("Después del fichero\n");

/*********************************************************/
/*********** Actualization (fluxes format) ***************/
/*********************************************************/
srand(1234);
if(schemeOrder == 1){
    while ((TMAX-t) > pow(10,-10)) {
        /*************************************************/
        /********************* FOU ***********************/
        /*************************************************/
        //-----------------------------------------------//
		/***************** CFL Condition *****************/
        if(equation == 2){
            uDT = -10000.0;
            for(i=0;i<NCELL;i++){
                if(fabs(u[i]) > uDT){
                    uDT = fabs(u[i]);
                }
            }
            DT = DX*CFL/uDT;
        }
        if(t + DT > TMAX){
            DT = TMAX - t + pow(10,-10);
        }
        /************* Spatial reconstruction ************/
		#ifdef openmp
			#pragma omp parallel for default(none) shared(u,u_R,u_L,NCELL)
		#endif
        for(i=0;i<NCELL-1;i++){
            u_L[i] = u[i];
            u_R[i] = u[i+1];
        }
        u_L[NCELL-1] = u[NCELL-1];
        u_R[NCELL-1] = u[0];
        /****************** Flux Calculus ****************/
        fluxCalculus(u, u_R, u_L, f_R, f_L, NCELL, lambda, equation);
        /******************* Source Term *****************/
        point_array = (int)(N_random*(double)(rand()/(RAND_MAX+1.0)));
        if(sourceTerm == 1){
			#ifdef openmp
				#pragma omp parallel for default(none) private(x1,x1_source,x2,x2_source,suma_source,n_source,index_rand) shared(N_random,point_array,xcell,DX,DT,amp_source,source,Nc,NCELL,L_source,randomarray)
            #endif
            for(i=0;i<NCELL;i++){
                x1          = xcell[i] + DX/2.0;
				x1_source   = 2.0*M_PI*x1/L_source;
                x2          = xcell[i] - DX/2.0;
				x2_source   = 2.0*M_PI*x2/L_source;
                suma_source = 0.0;
				n_source 	= 1;
                while(n_source<=Nc){
					if(point_array+n_source>=N_random){
						index_rand = point_array + n_source - N_random;
					}else{
						index_rand = point_array + n_source;
					}
					suma_source = suma_source + randomarray[index_rand]/pow(M_PI*fabs(n_source),1.5)*(sin(n_source*x1_source) - sin(n_source*x2_source));
					n_source = n_source + 1;
                }
                source[i]   = amp_source/pow(DT,0.5)*suma_source*L_source;
            }
        }
        /************* Actualization equation ************/
        unew[0] = u[0] - DT/DX*(f_L[0] - f_R[NCELL-1]) + sourceTerm*DT/DX*source[0];
		#ifdef openmp
			#pragma omp parallel for default(none) shared(unew,sourceTerm,DT,source,NCELL,DX,u,f_R,f_L)
        #endif
		for(i=1;i<NCELL;i++){
            unew[i] = u[i] - DT/DX*(f_L[i] - f_R[i-1]) + sourceTerm*DT/DX*source[i];
        }
        /****************** Data dumping *****************/
        for(i=0;i<NCELL;i++){
            u[i] = unew[i];
        }
		/****************** Data writing *****************/
		if(t > tthreshold){
			taux = taux + DT;
			if(taux > tdumping){
				sprintf (filename, "SOLUTION1dX_%d.txt",ntwrite);
				SF = fopen(filename, "wb");
				for(i=0;i<NCELL;i++){
					fprintf(SF,"%.16f\t %.16f\n",xcell[i],u[i]);
				}
				fclose(SF);
				taux    = 0.0;
				ntwrite = ntwrite + 1;
			}
		}
        t = t + DT;
    }
}else{
    while ((TMAX-t) > pow(10,-10)) {
        /*************************************************/
        /********************* RK3 ***********************/
        /*************************************************/
		//-----------------------------------------------//
		/***************** CFL Condition *****************/
        if(equation == 2){
            uDT = -10000.0;
            for(i=0;i<NCELL;i++){
                if(fabs(u[i]) > uDT){
                    uDT = fabs(u[i]);
                }
				//printf("uDT = %f\n",uDT);
            }
            DT = DX*CFL/uDT;
        }
        if(t + DT > TMAX){
            DT = TMAX - t + pow(10,-10);
        }
		/************* Spatial reconstruction ************/
        for(i=0;i<NCELL;i++){
            unew[i] = u[i];
        }
        for(k=0;k<schemeOrder;k++){
            switch(spatialOrder){
            case 3:
                rightLeft3(unew, 0, 1, 2, NCELL-1, 0, 1, u_R, u_L, NCELL, spatialScheme);
				#ifdef openmp
					#pragma omp parallel for default(none) shared(unew,NCELL,u_R,u_L,spatialScheme)				
                #endif
				for(i=1;i<NCELL-2;i++){
                    rightLeft3(unew, i, i+1, i+2, i-1, i, i+1, u_R, u_L, NCELL, spatialScheme);
                }
                rightLeft3(unew, NCELL-2, NCELL-1, 0, NCELL-3, NCELL-2, NCELL-1, u_R, u_L, NCELL, spatialScheme);
                rightLeft3(unew, NCELL-1,       0, 1, NCELL-2, NCELL-1,       0, u_R, u_L, NCELL, spatialScheme);
                break;
            case 5:
                rightLeft5(unew, NCELL-1, 0, 1, 2, 3, NCELL-2, NCELL-1, 0, 1, 2, u_R, u_L, NCELL, spatialScheme);
                rightLeft5(unew,       0, 1, 2, 3, 4, NCELL-1,       0, 1, 2, 3, u_R, u_L, NCELL, spatialScheme);
				#ifdef openmp
					#pragma omp parallel for default(none) shared(unew,NCELL,u_R,u_L,spatialScheme)						
                #endif
				for(i=2;i<NCELL-3;i++){
                    rightLeft5(unew, i-1, i, i+1, i+2, i+3, i-2, i-1, i, i+1, i+2, u_R, u_L, NCELL, spatialScheme);
                }
                rightLeft5(unew, NCELL-4, NCELL-3, NCELL-2, NCELL-1, 0, NCELL-5, NCELL-4, NCELL-3, NCELL-2, NCELL-1, u_R, u_L, NCELL, spatialScheme);
                rightLeft5(unew, NCELL-3, NCELL-2, NCELL-1,       0, 1, NCELL-4, NCELL-3, NCELL-2, NCELL-1,       0, u_R, u_L, NCELL, spatialScheme);
                rightLeft5(unew, NCELL-2, NCELL-1,       0,       1, 2, NCELL-3, NCELL-2, NCELL-1,       0,       1, u_R, u_L, NCELL, spatialScheme);
                break;
			case 7:
				rightLeft7(unew, NCELL-2, NCELL-1, 0, 1, 2, 3, 4, NCELL-3, NCELL-2, NCELL-1, 0, 1, 2, 3, u_R, u_L, NCELL, spatialScheme);
				rightLeft7(unew, NCELL-1, 	    0, 1, 2, 3, 4, 5, NCELL-2, NCELL-1,       0, 1, 2, 3, 4, u_R, u_L, NCELL, spatialScheme);
				rightLeft7(unew,       0,       1, 2, 3, 4, 5, 6, NCELL-1,       0,       1, 2, 3, 4, 5, u_R, u_L, NCELL, spatialScheme);
				for(i=3;i<NCELL-4;i++){
					rightLeft7(unew, i-2, i-1, i, i+1, i+2, i+3, i+4, i-3, i-2, i-1, i, i+1, i+2, i+3, u_R, u_L, NCELL, spatialScheme);
				}
				rightLeft7(unew, NCELL-6, NCELL-5, NCELL-4, NCELL-3, NCELL-2, NCELL-1, 0, NCELL-7, NCELL-6, NCELL-5, NCELL-4, NCELL-3, NCELL-2, NCELL-1, u_R, u_L, NCELL, spatialScheme);
				rightLeft7(unew, NCELL-5, NCELL-4, NCELL-3, NCELL-2, NCELL-1,       0, 1, NCELL-6, NCELL-5, NCELL-4, NCELL-3, NCELL-2, NCELL-1,       0, u_R, u_L, NCELL, spatialScheme);
				rightLeft7(unew, NCELL-4, NCELL-3, NCELL-2, NCELL-1,       0,       1, 2, NCELL-5, NCELL-4, NCELL-3, NCELL-2, NCELL-1,       0,       1, u_R, u_L, NCELL, spatialScheme);
				rightLeft7(unew, NCELL-3, NCELL-2, NCELL-1,       0,       1,       2, 3, NCELL-4, NCELL-3, NCELL-2, NCELL-1,       0,       1,       2, u_R, u_L, NCELL, spatialScheme);
				break;
            }
            /****************** Flux Calculus ****************/
            fluxCalculus(unew, u_R, u_L, f_R, f_L, NCELL, lambda, equation);
            /************* Actualization equation ************/
            unew[0] = a_eq[k]*u[0] + b_eq[k]*unew[0] - c_eq[k]*DT/DX*(f_L[0] - f_R[NCELL-1]);
			#ifdef openmp
				#pragma omp parallel for default(none) shared(a_eq,k,u,b_eq,unew,c_eq,DT,DX,f_R,f_L,NCELL)
			#endif
			for(i=1;i<NCELL;i++){
                unew[i] = a_eq[k]*u[i] + b_eq[k]*unew[i] - c_eq[k]*DT/DX*(f_L[i] - f_R[i-1]);
            }
        }
        /******************* Source Term *****************/
		point_array = (int)(N_random*(double)(rand()/(RAND_MAX+1.0)));
		//printf("point_array = %d\n",point_array);
        if(sourceTerm == 1){
			#ifdef openmp
				#pragma omp parallel for default(none) private(x1,x1_source,x2,x2_source,suma_source,n_source,index_rand) shared(N_random,point_array,xcell,DX,DT,amp_source,source,Nc,NCELL,L_source,randomarray)
            #endif
			for(i=0;i<NCELL;i++){
                x1          = xcell[i] + DX/2.0;
				x1_source   = 2*M_PI*x1/L_source;
                x2          = xcell[i] - DX/2.0;
				x2_source   = 2*M_PI*x2/L_source;
                suma_source = 0.0;
				n_source 	= 1;
                while(n_source<=Nc){
					if(point_array+n_source>=N_random){
						index_rand = point_array + n_source - N_random;
					}else{
						index_rand = point_array + n_source;
					}
					suma_source = suma_source + randomarray[index_rand]/pow(M_PI*fabs(n_source),1.5)*(sin(n_source*x1_source) - sin(n_source*x2_source));
					n_source = n_source + 1;
                }
                source[i]   = amp_source/pow(DT,0.5)*suma_source*L_source;
            }
        }
        /************* Actualization equation ************/
		#ifdef openmp
			#pragma omp parallel for default(none) shared(unew,sourceTerm,DT,source,NCELL,DX)
        #endif
		for(i=0;i<NCELL;i++){
            unew[i] = unew[i] + sourceTerm*DT/DX*source[i];
        }
        /****************** Data dumping *****************/
        for(i=0;i<NCELL;i++){
            u[i] = unew[i];
        }
		/****************** Data writing *****************/
		if(t > tthreshold){
			taux = taux + DT;
			if(taux > tdumping){
				sprintf (filename, "SOLUTION1dX_%d.txt",ntwrite);
				SF = fopen(filename, "wb");
				for(i=0;i<NCELL;i++){
					fprintf(SF,"%.16f\t %.16f\n",xcell[i],u[i]);
				}
				fclose(SF);
				taux    = 0.0;
				ntwrite = ntwrite + 1;
			}
		}
		/*
		if(t > 29.0){
			sprintf (filename, "SOLUTION1dX_%d.txt",ntwrite);
			SF = fopen(filename, "wb");
			for(i=0;i<NCELL;i++){
				fprintf(SF,"%.16f\t %.16f\n",xcell[i],u[i]);
			}
			fclose(SF);
			ntwrite = ntwrite + 1;
		}
		*/
		t = t + DT;
		nt = nt + 1;
		if(nt%1000000 == 0)
			printf("nt = %d  - DT = %f - t = %f - uDT = %f\n",nt,DT,t,uDT);
    }
	printf("Número de pasos temporales: %d\n",nt);
	
}
}

/*******************************************************/
/********************** Functions **********************/
/*******************************************************/

void fluxCalculus(double *vector, double *in_R, double *in_L, double *out_R, double *out_L, int N, double *lambda, int equation){
    int i, contador2;
    contador2 = 0;
    double u_mean[N], u_minus[N], u_plus[N], du[N], lambda_plus[N], lambda_minus[N], lambda_R[N], lambda_L[N];
    double lambda_mean[N];
    for(i=0;i<N;i++){
    u_mean[i]   	= 0.0;
    u_minus[i]  	= 0.0;
    u_plus[i] 		= 0.0;
	du[i] 	    	= 0.0;
	lambda_minus[i] = 0.0;
	lambda_plus[i]  = 0.0;
	lambda_R[i]     = 0.0;
	lambda_L[i]     = 0.0;
	lambda_mean[i]  = 0.0;
    }
    for(i=0;i<N-1;i++){ //by walls
        lambda_L[i] = lambda[i];
        lambda_R[i] = lambda[i+1];
    }
    lambda_L[N-1] = lambda[N-1];
    lambda_R[N-1] = lambda[0];
    for(i=0;i<N;i++){ // by walls
        du[i] = in_R[i] - in_L[i];
    }
    /************* Fluxes initialization ************/
    switch(equation){
    case 1:
        for(i=0;i<N;i++){
            out_R[i] = lambda[i]*in_R[i];
            out_L[i] = lambda[i]*in_L[i];
        }
        break;
    case 2:
        for(i=0;i<N;i++){
            out_R[i] = in_R[i]*in_R[i]/2.0;
            out_L[i] = in_L[i]*in_L[i]/2.0;

        }
        break;
    }
    /************* Fluxes actualization ************/
	switch(equation){
	case 1:
		for(i=0;i<N;i++){
            lambda_mean[i]  = (lambda_R[i] + lambda_L[i])/2.0;
			lambda_minus[i] = (lambda_mean[i] - fabs(lambda_mean[i]))/2.0;
			lambda_plus[i]  = (lambda_mean[i] + fabs(lambda_mean[i]))/2.0;
		}
		for(i=0;i<N;i++){
			if(lambda[i] > 0.0){
				out_R[i] = out_R[i] - lambda_plus[i]*du[i];
			} else {
				out_L[i] = out_L[i] + lambda_minus[i]*du[i];
			}
		}
		break;
	case 2:
		for(i=0;i<N;i++){
            u_mean[i]  = (in_R[i] + in_L[i])/2.0;
            // printf("u_mean[%d] = %f\n",i,u_mean[i]);
            if(in_L[i] < 0.0 && in_R[i] > 0.0){
                u_minus[i] = in_L[i]*(in_R[i] - u_mean[i])/(in_R[i] - in_L[i]);
                out_L[i]   = out_L[i] - u_minus[i]*du[i];
                u_plus[i]  = in_R[i]*(u_mean[i] - in_L[i])/(in_R[i] - in_L[i]);
                out_R[i]   = out_R[i] - u_plus[i]*du[i];
            }else{
                u_plus[i]  = (u_mean[i] + fabs(u_mean[i]))/2.0;
                u_minus[i] = (u_mean[i] - fabs(u_mean[i]))/2.0;
                if(u_mean[i] > 0.0){
                    out_R[i] = out_R[i] - u_plus[i]*du[i];
                }else{
                    out_L[i] = out_L[i] + u_minus[i]*du[i];
                }
            }
		}
        break;
	}
}

void rightLeft3(double *vector, int R_L, int R_C, int R_R, int L_L, int L_C, int L_R, double *out_R, double *out_L, int N, int spatialScheme){
    double g0[2], g1[2], g2[2], beta0[2], beta1[2], beta2[2], alpha0[2], alpha1[2], alpha2[2], w0[2], w1[2], w2[2], epsilon;
    int i;
    // Initialization
    for(i=0;i<2;i++){
        g0[i]     = 0.0;
        g1[i]     = 0.0;
        g2[i]     = 0.0;
        alpha0[i] = 0.0;
        alpha1[i] = 0.0;
        alpha2[i] = 0.0;
        w0[i]     = 0.0;
        w1[i]     = 0.0;
        w2[i]     = 0.0;
        beta0[i]  = 0.0;
        beta1[i]  = 0.0;
        beta2[i]  = 0.0;
    }
    epsilon = 0.0;
    epsilon = pow(10,-10);

    // Initial values
    g0[0] = 2.0/3.0;
    g0[1] = 1.0 - g0[0];
    g1[0] = 1.0/3.0;
    g1[1] = 1.0 - g1[0];

    switch(spatialScheme){
    case 1:
        alpha0[0] = g0[0];
        alpha0[1] = g0[1];
        alpha1[0] = g1[0];
        alpha1[1] = g1[1];
        break;
    case 2:
        beta0[0] = (vector[R_C]-vector[R_L])*(vector[R_C]-vector[R_L]);
        beta0[1] = (vector[L_C]-vector[L_L])*(vector[L_C]-vector[L_L]);
        beta1[0] = (vector[R_R]-vector[R_C])*(vector[R_R]-vector[R_C]);
        beta1[1] = (vector[L_R]-vector[L_C])*(vector[L_R]-vector[L_C]);
        alpha0[0] = g0[0]/((beta0[0]+epsilon)*(beta0[0]+epsilon));
        alpha0[1] = g0[1]/((beta0[1]+epsilon)*(beta0[1]+epsilon));
        alpha1[0] = g1[0]/((beta1[0]+epsilon)*(beta1[0]+epsilon));
        alpha1[1] = g1[1]/((beta1[1]+epsilon)*(beta1[1]+epsilon));
        break;
    }
    w0[0] = alpha0[0]/(alpha0[0]+alpha1[0]);
    w0[1] = alpha0[1]/(alpha0[1]+alpha1[1]);
    w1[0] = alpha1[0]/(alpha0[0]+alpha1[0]);
    w1[1] = alpha1[1]/(alpha0[1]+alpha1[1]);
    out_R[R_L] = w0[0]*(vector[R_C]/2.0  + vector[R_L]/2.0)     + w1[0]*(-vector[R_R]/2.0 + 3.0*vector[R_C]/2.0);
    out_L[L_C] = w0[1]*(-vector[L_L]/2.0 + 3.0*vector[L_C]/2.0) + w1[1]*(vector[L_C]/2.0  + vector[L_R]/2.0);
}

void rightLeft5(double *vector, int R_LL, int R_L, int R_C, int R_R, int R_RR, int L_LL, int L_L, int L_C, int L_R, int L_RR, double *out_R, double *out_L, int N, int spatialScheme){
    double g0[2], g1[2], g2[2], beta0[2], beta1[2], beta2[2], alpha0[2], alpha1[2], alpha2[2], w0[2], w1[2], w2[2], epsilon;
    int i;
    // Initialization
    for(i=0;i<2;i++){
        g0[i]     = 0.0;
        g1[i]     = 0.0;
        g2[i]     = 0.0;
        alpha0[i] = 0.0;
        alpha1[i] = 0.0;
        alpha2[i] = 0.0;
        w0[i]     = 0.0;
        w1[i]     = 0.0;
        w2[i]     = 0.0;
        beta0[i]  = 0.0;
        beta1[i]  = 0.0;
        beta2[i]  = 0.0;
    }
    epsilon = 0.0;
    // Initial values
    g0[0] = 3.0/10.0;
    g0[1] = 1.0/10.0;
    g1[0] = 3.0/5.0;
    g1[1] = 3.0/5.0;
    g2[0] = 1.0/10.0;
    g2[1] = 3.0/10.0;
    epsilon = pow(10,-14);
    switch(spatialScheme){
    case 1:
        alpha0[0] = g0[0];
        alpha0[1] = g0[1];
        alpha1[0] = g1[0];
        alpha1[1] = g1[1];
        alpha2[0] = g2[0];
        alpha2[1] = g2[1];
        break;
    case 2:
        beta0[0] = 13.0/12*(vector[R_LL] - 2*vector[R_L] + vector[R_C])*(vector[R_LL]  - 2*vector[R_L] + vector[R_C])  + 1.0/4*(vector[R_LL] - 4*vector[R_L] + 3*vector[R_C])*(vector[R_LL] - 4*vector[R_L] + 3*vector[R_C]);
        beta0[1] = 13.0/12*(vector[L_LL] - 2*vector[L_L] + vector[L_C])*(vector[L_LL]  - 2*vector[L_L] + vector[L_C])  + 1.0/4*(vector[L_LL] - 4*vector[L_L] + 3*vector[L_C])*(vector[L_LL] - 4*vector[L_L] + 3*vector[L_C]);
        beta1[0] = 13.0/12*(vector[R_L]  - 2*vector[R_C] + vector[R_R])*(vector[R_L]   - 2*vector[R_C] + vector[R_R])  + 1.0/4*(vector[R_L] - vector[R_R])*(vector[R_L] - vector[R_R]);
        beta1[1] = 13.0/12*(vector[L_L]  - 2*vector[L_C] + vector[L_R])*(vector[L_L]   - 2*vector[L_C] + vector[L_R])  + 1.0/4*(vector[L_L] - vector[L_R])*(vector[L_L] - vector[L_R]);
        beta2[0] = 13.0/12*(vector[R_C]  - 2*vector[R_R] + vector[R_RR])*(vector[R_C]  - 2*vector[R_R] + vector[R_RR]) + 1.0/4*(3*vector[R_C] - 4*vector[R_R] + vector[R_RR])*(3*vector[R_C] - 4*vector[R_R] + vector[R_RR]);
        beta2[1] = 13.0/12*(vector[L_C]  - 2*vector[L_R] + vector[L_RR])*(vector[L_C]  - 2*vector[L_R] + vector[L_RR]) + 1.0/4*(3*vector[L_C] - 4*vector[L_R] + vector[L_RR])*(3*vector[L_C] - 4*vector[L_R] + vector[L_RR]);
        alpha0[0] = g0[0]/((beta0[0] + epsilon)*(beta0[0] + epsilon));
        alpha0[1] = g0[1]/((beta0[1] + epsilon)*(beta0[1] + epsilon));
        alpha1[0] = g1[0]/((beta1[0] + epsilon)*(beta1[0] + epsilon));
        alpha1[1] = g1[1]/((beta1[1] + epsilon)*(beta1[1] + epsilon));
        alpha2[0] = g2[0]/((beta2[0] + epsilon)*(beta2[0] + epsilon));
        alpha2[1] = g2[1]/((beta2[1] + epsilon)*(beta2[1] + epsilon));
        break;
    }
    w0[0] = alpha0[0]/(alpha0[0] + alpha1[0] + alpha2[0]);
    w0[1] = alpha0[1]/(alpha0[1] + alpha1[1] + alpha2[1]);
    w1[0] = alpha1[0]/(alpha0[0] + alpha1[0] + alpha2[0]);
    w1[1] = alpha1[1]/(alpha0[1] + alpha1[1] + alpha2[1]);
    w2[0] = alpha2[0]/(alpha0[0] + alpha1[0] + alpha2[0]);
    w2[1] = alpha2[1]/(alpha0[1] + alpha1[1] + alpha2[1]);
    out_R[R_L] = w0[0]*(1/3.0*vector[R_C]  + 5.0/6*vector[R_L] - 1/6.0*vector[R_LL]) + w1[0]*(-1/6.0*vector[R_R] + 5/6.0*vector[R_C] + 1/3.0*vector[R_L]) + w2[0]*(1/3.0*vector[R_RR] - 7/6.0*vector[R_R] + 11/6.0*vector[R_C]);
    out_L[L_C] = w0[1]*(1/3.0*vector[L_LL] - 7/6.0*vector[L_L] + 11/6.0*vector[L_C]) + w1[1]*(-1/6.0*vector[L_L] + 5/6.0*vector[L_C] + 1/3.0*vector[L_R]) + w2[1]*(1/3.0*vector[L_C]  + 5/6.0*vector[L_R] - 1/6.0*vector[L_RR]);
}

void rightLeft7(double *vector, int R_LLL, int R_LL, int R_L, int R_C, int R_R, int R_RR, int R_RRR, int L_LLL, int L_LL, int L_L, int L_C, int L_R, int L_RR, int L_RRR, double *out_R, double *out_L, int N, int spatialScheme){
	double g0[2], g1[2], g2[2], g3[2], beta0[2], beta1[2], beta2[2], beta3[2], alpha0[2], alpha1[2], alpha2[2], alpha3[2], w0[2], w1[2], w2[2], w3[2], epsilon;
	int i;
    // Initialization
    for(i=0;i<2;i++){
        g0[i]     = 0.0;
        g1[i]     = 0.0;
        g2[i]     = 0.0;
		g3[i] 	  = 0.0;
		alpha0[i] = 0.0;
        alpha1[i] = 0.0;
        alpha2[i] = 0.0;
		alpha3[i] = 0.0;
        w0[i]     = 0.0;
        w1[i]     = 0.0;
        w2[i]     = 0.0;
		w3[i]     = 0.0;
        beta0[i]  = 0.0;
        beta1[i]  = 0.0;
        beta2[i]  = 0.0;
		beta3[i]  = 0.0;
    }
	epsilon = 0.0;
	// Initial values
	g0[0] = 4.0/35.0;
    g0[1] = 1.0/35.0;
	g1[0] = 18.0/35.0;
    g1[1] = 12.0/35.0;
    g2[0] = 12.0/35.0;
    g2[1] = 18.0/35.0;
    g3[0] = 1.0/35.0;
	g3[1] = 4.0/35.0;
	epsilon = pow(10,-14);
	switch(spatialScheme){
    case 1:
        alpha0[0] = g0[0];
        alpha0[1] = g0[1];
        alpha1[0] = g1[0];
        alpha1[1] = g1[1];
        alpha2[0] = g2[0];
        alpha2[1] = g2[1];
		alpha3[0] = g3[0];
        alpha3[1] = g3[1];
        break;
    case 2:
        beta0[0] = vector[R_LLL]*(547.0*vector[R_LLL] - 3882.0*vector[R_LL] + 4642.0*vector[R_L] - 1854.0*vector[R_C]) + vector[R_LL]*(7043.0*vector[R_LL] - 17246.0*vector[R_L] + 7042.0*vector[R_C]) + vector[R_L]*(11003.0*vector[R_L] - 9402.0*vector[R_C]) + vector[R_C]*2107.0*vector[R_C];
		beta0[1] = vector[L_LLL]*(547.0*vector[L_LLL] - 3882.0*vector[L_LL] + 4642.0*vector[L_L] - 1854.0*vector[L_C]) + vector[L_LL]*(7043.0*vector[L_LL] - 17246.0*vector[L_L] + 7042.0*vector[L_C]) + vector[L_L]*(11003.0*vector[L_L] - 9402.0*vector[L_C]) + vector[L_C]*2107.0*vector[L_C];
		beta1[0] = vector[R_LL]*(267.0*vector[R_LL] - 1642.0*vector[R_L] + 1602.0*vector[R_C] - 494.0*vector[R_R]) + vector[R_L]*(2843.0*vector[R_L] - 5966.0*vector[R_C] + 1922.0*vector[R_R]) + vector[R_C]*(3443.0*vector[R_C] - 2522.0*vector[R_R]) + vector[R_R]*547.0*vector[R_R];
		beta1[1] = vector[L_LL]*(267.0*vector[L_LL] - 1642.0*vector[L_L] + 1602.0*vector[L_C] - 494.0*vector[L_R]) + vector[L_L]*(2843.0*vector[L_L] - 5966.0*vector[L_C] + 1922.0*vector[L_R]) + vector[L_C]*(3443.0*vector[L_C] - 2522.0*vector[L_R]) + vector[L_R]*547.0*vector[L_R];
		beta2[0] = vector[R_L]*(547.0*vector[R_L] - 2522.0*vector[R_C] + 1922.0*vector[R_R] - 494.0*vector[R_RR]) + vector[R_C]*(3443.0*vector[R_C] - 5966.0*vector[R_R] + 1602*vector[R_RR]) + vector[R_R]*(2843.0*vector[R_R] - 1642*vector[R_RR]) + vector[R_RR]*267.0*vector[R_RR];
		beta2[1] = vector[L_L]*(547.0*vector[L_L] - 2522.0*vector[L_C] + 1922.0*vector[L_R] - 494.0*vector[L_RR]) + vector[L_C]*(3443.0*vector[L_C] - 5966.0*vector[L_R] + 1602*vector[L_RR]) + vector[L_R]*(2843.0*vector[L_R] - 1642*vector[L_RR]) + vector[L_RR]*267.0*vector[L_RR];
		beta3[0] = vector[R_C]*(2107.0*vector[R_C] - 9402.0*vector[R_R] + 7042.0*vector[R_RR] - 1854.0*vector[R_RRR]) + vector[R_R]*(11003.0*vector[R_R] - 17246.0*vector[R_RR] + 4642.0*vector[R_RRR]) + vector[R_RR]*(7043.0*vector[R_RR] - 3882.0*vector[R_RRR]) + vector[R_RRR]*547.0*vector[R_RRR];
		beta3[1] =vector[L_C]*(2107.0*vector[L_C] - 9402.0*vector[L_R] + 7042.0*vector[L_RR] - 1854.0*vector[L_RRR]) + vector[L_R]*(11003.0*vector[L_R] - 17246.0*vector[L_RR] + 4642.0*vector[L_RRR]) + vector[L_RR]*(7043.0*vector[L_RR] - 3882.0*vector[L_RRR]) + vector[L_RRR]*547.0*vector[L_RRR];
		alpha0[0] = g0[0]/((beta0[0] + epsilon)*(beta0[0] + epsilon));
        alpha0[1] = g0[1]/((beta0[1] + epsilon)*(beta0[1] + epsilon));
        alpha1[0] = g1[0]/((beta1[0] + epsilon)*(beta1[0] + epsilon));
        alpha1[1] = g1[1]/((beta1[1] + epsilon)*(beta1[1] + epsilon));
        alpha2[0] = g2[0]/((beta2[0] + epsilon)*(beta2[0] + epsilon));
        alpha2[1] = g2[1]/((beta2[1] + epsilon)*(beta2[1] + epsilon));
        alpha3[0] = g3[0]/((beta3[0] + epsilon)*(beta3[0] + epsilon));
        alpha3[1] = g3[1]/((beta3[1] + epsilon)*(beta3[1] + epsilon));
        break;
    }
	w0[0] = alpha0[0]/(alpha0[0] + alpha1[0] + alpha2[0] + alpha3[0]);
    w0[1] = alpha0[1]/(alpha0[1] + alpha1[1] + alpha2[1] + alpha3[1]);
    w1[0] = alpha1[0]/(alpha0[0] + alpha1[0] + alpha2[0] + alpha3[0]);
    w1[1] = alpha1[1]/(alpha0[1] + alpha1[1] + alpha2[1] + alpha3[1]);
    w2[0] = alpha2[0]/(alpha0[0] + alpha1[0] + alpha2[0] + alpha3[0]);
    w2[1] = alpha2[1]/(alpha0[1] + alpha1[1] + alpha2[1] + alpha3[1]);
	w3[0] = alpha3[0]/(alpha0[0] + alpha1[0] + alpha2[0] + alpha3[0]);
    w3[1] = alpha3[1]/(alpha0[1] + alpha1[1] + alpha2[1] + alpha3[1]);
	out_R[R_L] = w0[0]*(1.0/4.0*vector[R_C]  + 13.0/12.0*vector[R_L] - 5.0/12.0*vector[R_LL] + 1.0/12.0*vector[R_LLL]) + w1[0]*(-1.0/12.0*vector[R_R] + 7.0/12.0*vector[R_C] + 7.0/12.0*vector[R_L] - 1.0/12.0*vector[R_LL]) + w2[0]*(1.0/12.0*vector[R_RR] - 5.0/12.0*vector[R_R] + 13.0/12.0*vector[R_C] + 1.0/4.0*vector[R_L]) + w3[0]*(-1.0/4.0*vector[R_RRR] + 13.0/12.0*vector[R_RR] - 23.0/12.0*vector[R_R] + 25.0/12.0*vector[R_C]);
    out_L[L_C] = w0[1]*(-1.0/4.0*vector[L_LLL] + 13.0/12.0*vector[L_LL] - 23.0/12.0*vector[L_L] + 25.0/12.0*vector[L_C]) + w1[1]*(1.0/12.0*vector[L_LL] - 5.0/12.0*vector[L_L] + 13.0/12.0*vector[L_C] + 1.0/4.0*vector[L_R]) + w2[1]*(-1.0/12.0*vector[L_L] + 7.0/12.0*vector[L_C]  + 7.0/12.0*vector[L_R] - 1.0/12.0*vector[L_RR]) + w3[1]*(1.0/4.0*vector[L_C] + 13.0/12.0*vector[L_R] - 5.0/12.0*vector[L_RR] + 1.0/12.0*vector[L_RRR]);
}
