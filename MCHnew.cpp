//author: Wilson Wang, 09/2016
#include <stdio.h>
#include <stdlib.h> /* for rand() */
//#include <unistd.h> /* for getpid() */
#include <time.h> /* for time() */
#include <math.h>
#include <assert.h>
#include<iostream>
#include <cstdlib>
#include <ctime>
#include <fstream>


//typedef struct {
//double x, y, z;
//} spin;

using namespace std;
//int nSteps = 1000;//int nSteps = 1000;
double J = 1.0, A = 1.0,D,a=1.0;
#define N0 999//rand()%(N0 + 1)/(float)(N0 + 1)
// (rand()%(N0 + 1)/(float)(N0 + 1)-0.5)
#define pi 3.1415926
#define ABS(a) ((a)<0.0?(-(a)):(a))
double rin3[20][20] = { 0 }, rin5[20][20] = { 0 }, rij[11][11] = { 0 },r_ij[20][20] = { 0 };

double spins[20][20][3] = { 0 };




//-------------------------------------relax begin------------------------------------------------------------------------------------
void relax(int Nt, int N, float D, float T) {
	int n, p, i, j, i0, j0, si, sj, sk, ip, im, jp, jm, ip0, im0, jp0, jm0;
	double  local01, local00;
	double s[3] = { 0 };
	double x, M = 0, beta = 1.0, cv, xi, E = 0, mu = 1.0, deltaE;
	double rc, ex = 0, exd = 0, exd1 = 0, ex01, exd01 = 0, exd101 = 0, rij, r_3, r_5, anisotropy0 = 0, anisotropy1 = 0, dipolar0 = 0, dipolar1 = 0;
	double  exo = 0, anisotropy0o = 0, rijo, r_3o, r_5o, exd1o = 0, exdo = 0, dipolar0o = 0, z, exo2 = 0, anisotropy0o2 = 0, exdo2 = 0, exd1o2 = 0, dipolar0o2 = 0;
	int ipo, imo, jpo, jmo, im0o, jm0o;
	rc = 5.0;
	double *theta = new double[N*N];
	double *phi = new double[N*N];
	//---------------------------------------------------------------------------------------------------------------------
	//begin ------------------------------------------------------------------------------------------------------
	for (n = 0; n<Nt; n++)
	{
		si = (int)(N*(rand() / (1.0*RAND_MAX))); //choose the random point 
		sj = (int)(N*(rand() / (1.0*RAND_MAX)));
		if (si >= N) si = N - 1; if (sj >= N) sj = N - 1;
		s[0] = spins[si][sj][0]; s[1] = spins[si][sj][1]; s[2] = spins[si][sj][2];
		i0 = si;
		ip = i0 + 1; if (ip >= N) ip = 0;
		im = i0 - 1; if (im<0) im = N - 1;
		j0 = sj;
		jp = j0 + 1; if (jp >= N) jp = 0;
		jm = j0 - 1; if (jm<0) jm = N - 1;
		//caculate the local energy----------------
		ex = (-s[0] * J*(spins[ip][sj][0] + spins[im][sj][0] + spins[si][jp][0] + spins[si][jm][0])) + (-s[1] * J*(spins[ip][sj][1] + spins[im][sj][1] + spins[si][jp][1] + spins[si][jm][1])) + (-s[2] * J*(spins[ip][sj][2] + spins[im][sj][2] + spins[si][jp][2] + spins[si][jm][2]));
		anisotropy0 = A*spins[si][sj][2] * spins[si][sj][2];//caculate the anisotropy
		im0 = i0 - 5; if (im0<0) im0 += N;//caculate the dipolar 
		jm0 = j0 - 5; if (jm0<0) jm0 += N;
		for (int i = 0; i <= 10; i++){
			im0 = im0 + 1; if (im0 >= N) im0 = 0;
			for (int j = 0; j <= 10; j++)
			{
				jm0 = jm0 + 1; if (jm0 >= N) jm0 = 0;
				rij = sqrt((jm0 - j0)*(jm0 - j0) + (im0 - i0)*(im0 - i0));

				if (rij <= rc&rij>0){
					r_3 = pow(rij, -3);
					r_5 = pow(rij, -5);
					exd = exd + r_3*(s[0] * (spins[im0][jm0][0]) + s[1] * (spins[im0][jm0][1]) + s[2] * spins[im0][jm0][2]);
					exd1 = exd1 + (-3 * r_5*((im0 - i0)*spins[si][sj][0] + (jm0 - j0)*spins[si][sj][1])*((im0 - i0)*spins[im0][jm0][0] + (jm0 - j0)*spins[im0][jm0][1]));//(s_i\cdot r_ij)(s_j\cdot r_ij)


				}
			}
		}
		dipolar0 = exd + exd1;
		local00 = 0.5*ex + anisotropy0 + 0.5*D*dipolar0;
		exd = 0, exd1 = 0;
		//angular change-----------------------------------------------------------------------------
		theta[si*N + sj] = 2.0*pi*(rand() / (1.0*RAND_MAX)); // % (N0 + 1) / (double)(N0 + 1));
		phi[si*N + sj] = 1.0*pi*(rand() / (1.0*RAND_MAX));
		s[0] = sin(theta[si*N + sj])*cos(phi[si*N + sj]);
		s[1] = sin(theta[si*N + sj])*sin(phi[si*N + sj]);
		s[2] = cos(theta[si*N + sj]);

		//s[1] = spins[si][sj][1]; s[2] = spins[si][sj][2]; s[3] = spins[si][sj][3];
		i0 = si;
		ip = i0 + 1; if (ip >= N) ip = 0;
		im = i0 - 1; if (im<0) im = N - 1;
		j0 = sj;
		jp = j0 + 1; if (jp >= N) jp = 0;
		jm = j0 - 1; if (jm<0) jm = N - 1;

		ex01 = (-s[0] * J*(spins[ip][sj][0] + spins[im][sj][0] + spins[si][jp][0] + spins[si][jm][0])) + (-s[1] * J*(spins[ip][sj][1] + spins[im][sj][1] + spins[si][jp][1] + spins[si][jm][1])) + (-s[2] * J*(spins[ip][sj][2] + spins[im][sj][2] + spins[si][jp][2] + spins[si][jm][2]));
		anisotropy1 = A*s[2] * s[2];
		im0 = i0 - 5; if (im0<0) im0 += N;
		jm0 = j0 - 5; if (jm0<0) jm0 += N;
		for (int i = 0; i <= 10; i++){
			im0 = im0 + 1; if (im0 >= N) im0 = 0;
			for (int j = 0; j <= 10; j++)
			{
				jm0 = jm0 + 1; if (jm0 >= N) jm0 = 0;
				//rij = sqrt(1.0*(j-5)*(j- 5) + 1.0*(i - 5)*(i - 5));
				rij = sqrt(1.0*(jm0 - j0)*(jm0 - j0) + 1.0*(im0 - i0)*(im0 - i0));
				if (rij <= rc&&rij>0){
					r_3 = pow(rij, -3);
					//r_3 = rin3[ABS(i - 5)][ABS(j - 5)];
					//r_5=rin5[ABS(i - 5)][ABS(j - 5)];
					r_5 = pow(rij, -5);
					exd01 = exd01 + r_3*(s[0] * (spins[im0][jm0][0]) + s[1] * (spins[im0][jm0][1]) + s[2] * spins[im0][jm0][2]);
					exd101 = exd101 + (-3 * r_5*((im0 - i0)*s[0] + (jm0 - j0)*s[1])*((im0 - i0)*spins[im0][jm0][0] + (jm0 - j0)*spins[im0][jm0][1]));//(s_i\cdot r_ij)(s_j\cdot r_ij)


				}
			}
		}
		dipolar1 = exd01 + exd101;
		local01 = 0.5*ex01 + anisotropy1 + 0.5*D*dipolar1;
		exd01 = 0; exd101 = 0;
		//printf("%f,%f\n", local01, local00);
		deltaE = local01 - local00;
		//printf("%f\n",deltaE);

		if (deltaE <0){
			spins[si][sj][0] = s[0];   spins[si][sj][1] = s[1];  spins[si][sj][2] = s[2];
		}
		else {

			x = (rand() / (1.0*RAND_MAX));
			if (exp(-deltaE / T) >x)
			{
				spins[si][sj][0] = s[0];   spins[si][sj][1] = s[1];  spins[si][sj][2] = s[2];
			}
		}

	}
	
	delete theta;
	delete phi;
}



//-------------------------------------------relax end--------------------------------------------------------------------

void MCsteps(int Nt,int N,float D,float T) {
	int n, p, i, j, i0, j0, si, sj, sk, ip, im, jp, jm, ip0, im0, jp0, jm0;
	double  local01, local00,avM, avM2, avE, avE2;
	double s[3] = {0};
	double x, M = 0, beta = 1.0, cv, xi, E = 0, mu = 1.0, deltaE;
	double rc, ex=0, exd = 0, exd1 = 0, ex01, exd01 = 0, exd101 = 0, rij, r_3, r_5, anisotropy0 = 0, anisotropy1 = 0, dipolar0 = 0, dipolar1 = 0;
	double  exo = 0, anisotropy0o = 0, rijo, r_3o, r_5o, exd1o = 0, exdo = 0, dipolar0o = 0, z, exo2 = 0, anisotropy0o2 = 0, exdo2 = 0, exd1o2 = 0, dipolar0o2 = 0;
	int ipo, imo, jpo, jmo, im0o, jm0o;
	
	rc = 5.0;
	
	
	double *theta = new double[N*N];
	double *phi = new double[N*N];

//	int N = 20;
//	spin *spins0;
//	spins0 = new spin[N*N];

//	double sijx = spins0[i + N*j].x;


	


	//printf("E %f\n", E);
	//exo = 0; anisotropy0o = 0;

	avM = 1.0*M, avM2 = 1.0*M*M, avE = 1.0*E, avE2 = 1.0*E*E;
	ofstream outfile;
	outfile.open("results.txt", ios::out | ios::app);
	if (!outfile){
		cout << "file cannot open.\n";
		abort();
	}

	//---------------------------------------------------------------------------------------------------------------------
	//begin the MC steps-------------------------------------------------------------------------------------------------------
	for (n = 0; n<Nt; n++) 
	{
		si = (int)(N*(rand()/(1.0*RAND_MAX))); //choose the random point 
		sj = (int)(N*(rand()/(1.0*RAND_MAX)));
		if (si >=N) si = N - 1; if (sj >=N) sj = N - 1;
		s[0] = spins[si][sj][0]; s[1] = spins[si][sj][1]; s[2] = spins[si][sj][2];
		i0 = si;
		ip = i0 + 1; if (ip >=N) ip = 0;
		im = i0 - 1; if (im<0) im = N - 1;
		j0 = sj;
		jp = j0 + 1; if (jp >= N) jp = 0;
		jm = j0 - 1; if (jm<0) jm = N - 1;
		//caculate the local energy----------------
		ex = (-s[0] * J*(spins[ip][sj][0] + spins[im][sj][0] + spins[si][jp][0] + spins[si][jm][0])) + (-s[1] * J*(spins[ip][sj][1] + spins[im][sj][1] + spins[si][jp][1] + spins[si][jm][1])) + (-s[2] * J*(spins[ip][sj][2] + spins[im][sj][2] + spins[si][jp][2] + spins[si][jm][2]));
		anisotropy0 = A*spins[si][sj][2] * spins[si][sj][2];//caculate the anisotropy
		im0 = i0 - 5; if (im0<0) im0 += N ;//caculate the dipolar 
		jm0 = j0 - 5; if (jm0<0) jm0 += N ;
		for (int i = 0; i <= 10; i++){
			im0 = im0 +1; if (im0>=N) im0 =0;
			for (int j = 0; j <= 10; j++)
			{
				jm0 = jm0 +1; if (jm0>=N) jm0 =0;
				rij = sqrt((jm0 - j0)*(jm0 - j0) + (im0 - i0)*(im0 - i0));
				
				if (rij <= rc&rij>0){
					r_3 = pow(rij, -3);
					r_5 = pow(rij, -5);
					exd=exd+r_3*(s[0] *(spins[im0][jm0][0]) +s[1] *(spins[im0][jm0][1]) +s[2] *spins[im0][jm0][2]);
					exd1 = exd1 + (-3 * r_5*((im0 - i0)*spins[si][sj][0] + (jm0 - j0)*spins[si][sj][1])*((im0 - i0)*spins[im0][jm0][0] + (jm0 - j0)*spins[im0][jm0][1]));//(s_i\cdot r_ij)(s_j\cdot r_ij)
					

				}
			}
		}
		dipolar0 = exd + exd1;
		local00 = 0.5*ex + anisotropy0 + 0.5*D*dipolar0;
		exd = 0, exd1 = 0;
		//angular change-----------------------------------------------------------------------------
		z = 1;// (rand() % (N0 + 1) / (double)(N0 + 1) - 0.5 > 0) ? 1 : -1;
		theta[si*N + sj] = 2.0*pi*(rand() / (1.0*RAND_MAX)); // % (N0 + 1) / (double)(N0 + 1));
		phi[si*N + sj] = 1.0*pi*(rand() / (1.0*RAND_MAX));
		s[0] = z*sin(theta[si*N + sj])*cos(phi[si*N + sj]);
		s[1] = z*sin(theta[si*N + sj])*sin(phi[si*N + sj]);
		s[2] = z*cos(theta[si*N + sj]);

		//s[1] = spins[si][sj][1]; s[2] = spins[si][sj][2]; s[3] = spins[si][sj][3];
		i0 = si;
		ip = i0 + 1; if (ip >= N) ip = 0;
		im = i0 - 1; if (im<0) im = N - 1;
		j0 = sj;
		jp = j0 + 1; if (jp >= N) jp = 0;
		jm = j0 - 1; if (jm<0) jm = N - 1;

		ex01 = (-s[0] * J*(spins[ip][sj][0] + spins[im][sj][0] + spins[si][jp][0] + spins[si][jm][0])) + (-s[1] * J*(spins[ip][sj][1] + spins[im][sj][1] + spins[si][jp][1] + spins[si][jm][1])) + (-s[2] * J*(spins[ip][sj][2] + spins[im][sj][2] + spins[si][jp][2] + spins[si][jm][2]));
		anisotropy1 = A*s[2] * s[2];



		/*distance
		for (j = -5; j <= 5; j++) {
			jnew0 = j0 + j; if (jnew0 >= N) jnew0 -= N; else if (jnew0 < 0) jnew0 += N;
			for (i = -5; i <= 5; i++) {
				inew0 = i0 + i; if (inew0 >= N) inew0 -= N; else if (inew0 < 0) inew0 += N;

				rij[i][j] = sqrt(j*j+i*i);

			}
		}


		*/





		im0 = i0 - 5; if (im0<0) im0 += N ;
		jm0 = j0 - 5; if (jm0<0) jm0 += N;
		for (int i = 0; i <= 10; i++){
			im0 = im0 + 1; if (im0 >= N) im0 = 0;
			for (int j = 0; j <= 10; j++)
			{
				jm0 = jm0 + 1; if (jm0 >= N) jm0 = 0;
				rij = sqrt(1.0*(j-5)*(j- 5) + 1.0*(i - 5)*(i - 5));
				//rij = sqrt(1.0*(jm0 - j0)*(jm0 - j0) + 1.0*(im0 - i0)*(im0 - i0));
				if (rij <= rc&&rij>0){
					r_3 = pow(rij, -3);
					//r_3 = rin3[ABS(i - 5)][ABS(j - 5)];
					//r_5=rin5[ABS(i - 5)][ABS(j - 5)];
					r_5 = pow(rij, -5);
					exd01 = exd01 + r_3*(s[0] * (spins[im0][jm0][0]) +s[1] * (spins[im0][jm0][1]) +s[2] * spins[im0][jm0][2]);
					exd101 = exd101 + (-3 * r_5*((im0 - i0)*s[0] + (jm0 - j0)*s[1])*((im0 - i0)*spins[im0][jm0][0] + (jm0 - j0)*spins[im0][jm0][1]));//(s_i\cdot r_ij)(s_j\cdot r_ij)


				}
			}
		}
		dipolar1 = exd01 + exd101;
		local01 = 0.5*ex01 + anisotropy1 + 0.5*D*dipolar1;
		exd01 = 0; exd101 = 0;
		//printf("%f,%f\n", local01, local00);
		deltaE = local01 - local00;
		//printf("%f\n",deltaE);

		if (deltaE < 0){ 
			spins[si][sj][0] = s[0];   spins[si][sj][1] = s[1];  spins[si][sj][2] = s[2]; 
		}
		else {

			x = (rand()/(1.0*RAND_MAX));
			if (exp(-deltaE/T) >x)
			{
				spins[si][sj][0] = s[0];   spins[si][sj][1] = s[1];  spins[si][sj][2] = s[2];
			}
		}


		//-------------------------------calculate energy and magnetism ----------------------------------------------------------------

		E = 0; M = 0; exdo2 = 0; exd1o2 = 0; anisotropy0o2 = 0; dipolar0o2 = 0;
		for (int i = 0; i < N; i++)
		{
			for (int j = 0; j < N; j++)
			{

				M = M + mu*spins[i][j][2];//magnetism


				ipo = i + 1; if (ipo >= N) ipo = 0;
				imo = i - 1; if (imo<0) imo = N - 1;
				jpo = j + 1; if (jpo >= N) jpo = 0;
				jmo = j - 1; if (jmo<0) jmo = N - 1;
				//caculate the  energy----------------
				exo2 = exo2+(-spins[i][j][0] * J*(spins[ipo][j][0] + spins[imo][j][0] + spins[i][jpo][0] + spins[i][jmo][0])) + (-spins[i][j][1] * J*(spins[ipo][j][1] + spins[imo][j][1] + spins[i][jpo][1] + spins[i][jmo][1])) + (-spins[i][j][2] * J*(spins[ipo][j][2] + spins[imo][j][2] + spins[i][jpo][2] + spins[i][jmo][2]));
				//printf("exo %f\n", exo);exo2=0,anisotropy0o2=0,exdo2=0,exd1o2=0,dipolar0o2=0;
				anisotropy0o2 = anisotropy0o2+A*spins[i][j][2] * spins[i][j][2];//caculate the anisotropy
				//printf("anisotropy %f\n", anisotropy0o);
				im0o = i - 5; if (im0o<0) im0o += N;//caculate the dipolar 
				jm0o = j - 5; if (jm0o<0) jm0o += N;
				for (int io = 0; io <= 10; io++)
				{
					im0o = im0o + 1; if (im0o >= N) im0o = 0;
					for (int jo = 0; jo <= 10; jo++)
					{
						jm0o = jm0o + 1; if (jm0o >= N) jm0o = 0;
						rijo = sqrt((jm0o - j)*(jm0o - j) + (im0o - i)*(im0o - i));
						//printf("rijo %lf\n", rijo);
						//printf("jm0o %d\n", im0o);
						//printf("i j %d,%d\n", i,j);


						if (rijo <= rc&&rijo>0){
							r_3o = pow(rijo, -3);
							r_5o = pow(rijo, -5);
							//printf("rc-5 %lf\n", pow(rc,-5));
							//printf("r-5 %lf\n", r_5o);
							//printf("exdo %lf\n", exdo);
							exdo2 = exdo2 + r_3o*(spins[i][j][0] * (spins[im0o][jm0o][0]) + spins[i][j][1] * (spins[im0o][jm0o][1]) + spins[i][j][2] * spins[im0o][jm0o][2]);
							//printf("exdo2 %lf\n", exdo2);
							exd1o2 = exd1o2 + (-3 * r_5o*((im0o - i)*spins[i][j][0] + (jm0o - j)*spins[i][j][1])*((im0o - i)*spins[im0o][jm0o][0] + (jm0o - j)*spins[im0o][jm0o][1]));//(s_i\cdot r_ij)(s_j\cdot r_ij)
							//printf("exd1o2 %lf\n", exd1o2);

						}
					}
				}
				dipolar0o2 = exdo2 + exd1o2;
				exdo2 = 0; exd1o2 = 0;
				E = E + 0.5*exo2 + anisotropy0o2 + 0.5*D*dipolar0o2;
				dipolar0o2 = 0;
				//printf("E %lf\n", E);
			}
		}
		//-------------------------------calculate energy and magnetism -----end-----------------------
		
		avM += 1.0*M;
		avM2 += 1.0*M*M;
		avE2 += 1.0*E*E;
		avE += 1.0*E;
		//printf("E %f,avE %f\n", E,avE);
		//if(vis) {if(n%visfreq==0) redraw();}
		E = 0; M = 0;
	}

	//--------------------------------------------------------------------------------------------------------------------
	//------ the MC steps-----end--------------------------------------------------------------------------------------------------
	avM = avM / (1.0*(Nt + 1));
	avM2 = avM2 / (1.0*(Nt + 1));
	avE = avE / (1.0*(Nt + 1));
	avE2 = avE2 / (1.0*(Nt + 1));
	cv = beta*(avE2 - avE*avE) / T;
	E = 0; M = 0;
	xi = beta*(avM2-avM*avM);
	//printf("E %f,avE %f\n", E, avE);
	//printf("T %f,cv %f,xi %f,avM %f,avM2 %f,avE %f,avE2, %f \n",T,cv,xi,avM,avM2,avE,avE2);
	printf("T %f,cv %f,xi %f,avE %f \n", T, cv, xi, avE);
	
	//printf("%f,%f\n", avE, T);
	outfile << T << "\t" <<cv<< "\t" << xi <<"\t"<< avM << "\t " << avE << "\t " << "\n" << endl;
	outfile.close();
	delete theta;
	delete phi;
}



//-------------------------main()---------------------------------------------------
int main(int argc, char *argv[]){
	double J = 1, A = 1,T;
	double a = 1.0;
	double rc = 5 * a;
	int N = 20;
	srand(time(NULL));

	double *theta = new double[N*N];
	double *phi = new double[N*N];

	for (int i = 0; i < N; i++){//initialize the system
		for (int j = 0; j < N; j++)
		{
			theta[i*N + j] = 2.0*pi*(rand() / (1.0*RAND_MAX));
			phi[i*N + j] = pi*(rand() / (1.0*RAND_MAX));
			spins[i][j][0] = sin(theta[i*N + j])*cos(phi[i*N + j]);
			spins[i][j][1] = sin(theta[i*N + j])*sin(phi[i*N + j]);
			spins[i][j][2] = cos(theta[i*N + j]);
			//printf("%f,%f,%f\n", spins[i][j][0], spins[i][j][1], spins[i][j][2]);

		}
	}

	for (int i = 0; i < 11; i++) //distance
	{
		for (int i = 0; i < 11; i++)
		{
			rij[i][j] = sqrt(i*i + j*j);
			rin3[i][j] = pow(rij[i][j], -3);
			rin5[i][j] = pow(rij[i][j], -5);
		}
	
	
	
	
	}
	//printf("%f,%f,%f\n", spins[10][10][0], spins[10][10][1], spins[10][10][2]);

	

	T = 0.2;
	for (int i = 0; i < 10; i++){

		T = T + 0.01;
		relax(500, 20, 0.2, T);
		MCsteps(1000, 20, 0.2, T);//MCsteps(int Nt,int N,float D,float T)
		//MCsteps(6000, 20, 0.1, 1);
	}
	//printf("...........................................\n");

	



	delete theta;
	delete phi;
	system("pause");

	//return 0;
}
