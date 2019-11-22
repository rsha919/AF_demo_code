
#define _CRT_SECURE_NO_WARNINGS
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#define X   100			
#define Y   6			
#define Z   6			
#define D1A   30.0				
#define D2A   3.0			
#define DDA   27.0	
#define dx  0.03
double  Y0[X+1][Y+1][Z+1];
double  Y1[X+1][Y+1][Z+1];
double  Y2[X+1][Y+1][Z+1];
double  Y3[X+1][Y+1][Z+1];
double Y10[X+1][Y+1][Z+1];
double Y11[X+1][Y+1][Z+1];
double Y12[X+1][Y+1][Z+1];
double Y13[X+1][Y+1][Z+1];
double Y14[X+1][Y+1][Z+1];
double Y15[X+1][Y+1][Z+1];
double Y16[X+1][Y+1][Z+1];
double Y17[X+1][Y+1][Z+1];
double Y18[X+1][Y+1][Z+1];
double Y22[X+1][Y+1][Z+1];
double Y24[X+1][Y+1][Z+1];
double Y25[X+1][Y+1][Z+1];
double Y26[X+1][Y+1][Z+1];
double Y27[X+1][Y+1][Z+1];
double Y28[X+1][Y+1][Z+1];
double Y29[X+1][Y+1][Z+1];
double Y30[X+1][Y+1][Z+1];
double Y31[X+1][Y+1][Z+1];
double Y32[X+1][Y+1][Z+1];
double dY0;
double dY1;
double dY2;
double dY3;
double dY10;
double dY11;
double dY12;
double dY13;
double dY14;
double dY15;
double dY16;
double dY17;
double dY18;
double dY22;
double dY24;
double dY25;
double dY26;
double dY27;
double dY28;
double dY29;
double dY30;
double dY31;
double dY32;
double Cao;   // millimolar (in Ionic_values)
double Ki;   // millimolar (in Ionic_values)
double Ko;   // millimolar (in Ionic_values)
double Nao;   // millimolar (in Ionic_values)
double F;   // coulomb_per_mole (in Membrane)
double R2;   // joule_per_kilomole_kelvin (R in Membrane)
double T;   // kelvin (in Membrane)
double ACh;   // millimolar (in Rate_modulation_experiments)
double ACh_on;   // dimensionless (in i_KACh)
double E_Ca;   // millivolt (in Ionic_values)
double E_K;   // millivolt (in Ionic_values)
double E_Na;   // millivolt (in Ionic_values)
double RTonF;   // millivolt (in Membrane)
double V;   // millivolt (in Membrane)
double i_tot;   // nanoA (in Membrane)
double Nai;   // millimolar (in Nai_concentration)
double i_CaL;   // nanoA (in i_CaL)
double i_KACh;   // nanoA (in i_KACh)
double i_Kr;   // nanoA (in i_Kr)
double i_Ks;   // nanoA (in i_Ks)
double i_Kur;   // nanoA (in i_Kur)
double i_NaCa;   // nanoA (in i_NaCa)
double i_NaK;   // nanoA (in i_NaK)
double alpha_h;   // per_second (in i_Na_h_gate)
double beta_h;   // per_second (in i_Na_h_gate)
double h_infinity;   // dimensionless (in i_Na_h_gate)
double tau_h;   // second (in i_Na_h_gate)
double alpha_m;   // per_second (in i_Na_m_gate)
double beta_m;   // per_second (in i_Na_m_gate)
double m_infinity;   // dimensionless (in i_Na_m_gate)
double tau_m;   // second (in i_Na_m_gate)
double i_Na;   // nanoA (in i_Na)
double i_to;   // nanoA (in i_to)
int g[X+1][Y+1][Z+1];
int h[X+1][Y+1][Z+1];
double xx[X+1][Y+1][Z+1];
double yy[X+1][Y+1][Z+1];
double zz[X+1][Y+1][Z+1];
double dc[X+1][Y+1][Z+1][10];
double df[X+1][Y+1][Z+1][10];
double cmdnbar = 0.050;   
double kmnai= 10.0;    
double kmko = 1.5;    
double kmnancx = 87.5;  
double kmcancx = 1.38;  
double kmpca = 0.0005;
double tautr = 0.180; 
double vjsr  = 0.009648e-5;
double vnsr = 0.110952e-5;
double vmyo = 1.3668e-5;
double trpnbar = 0.070;   
double kmtrpn = 0.0005;  
double kmcmdn = 0.00238;  
double tau_u = 0.008;
double C_CRN = 0.0001; 
double gks = 0.0129;
double gito = 0.01652;
double gnab = 6.74e-5; 
double gcab = 1.13e-4; 
double gkb = 0.0; 
double gna = 0.78; 
double gcalbar = 0.01238; 
double gkr = 0.00294; 
double gk1 = 0.009; 
double gkurMultiplyer = 0.001; 
double gkach = 0.0135; 
double ibarpaca = 160.0;  
double ibarpca = 0.0275; 
double ibarnak = 0.06;
double cmdn;
double ua;
double ui;
double xs;
double a;
double igate;
double m;
double hgate;
double jgate;
double d;
double f;
double fca;
double xr;
double u;
double v;
double w;
double i_rel;
double i_up;
double i_leak;
double i_tr;
double Na_ion_tot;
double K_ion_tot;
double Ca_ion_tot;
double b1Cai;
double b2Cai;
double fn;
double u_infinity;
double tau_v;
double v_infinity;
double tau_w;
double w_infinity;
double alpha_ua;
double beta_ua;
double tau_ua;
double ua_infinity;
double alpha_ui;
double beta_ui;
double tau_ui;
double ui_infinity;
double alpha_xs;
double beta_xs;
double tau_xs;
double xs_infinity;
double alpha_a;
double beta_a;
double tau_a;
double a_infinity;
double alpha_i;
double beta_i;
double tau_i;
double i_infinity;
double alpha_j;
double beta_j;
double j_infinity;
double tau_j;
double d_infinity;
double tau_d;
double f_infinity;
double tau_f;
double fca_infinity;
double tau_fca;
double alpha_xr;
double beta_xr;
double tau_xr;
double xr_infinity;
double i_K1;
double i_pCa;
double i_bNa;
double i_bCa;
double i_bK;
double gkur;
double alpha_yACh;
double beta_yACh;
double tau_yACh;
double yACh_infinity;
double yACh;
double D1;
double D2;
double DD;

int main (int argc, char **argv)
{

	int x,y,z, i,j,k;
	int num;
	int isbound;
	int gg[27];
	float root2 = sqrt(2.0);
	float root3 = sqrt(3.0);
	float ic, ir, il, imax;
	float tflt;
	double du;
	double dudx, dudy, dudz;
	double dudx2, dudy2, dudz2;
	double dudxdy,dudxdz, dudydz;
	char *str;
	char *str2;
	char *str3;
	FILE *out;
	FILE *in;
	double t, udt;    
	double end_time;
	double steps; 
	long int increment;
	int incrementDivisor = 400;
	double ptstim = 0.0;
    double counter = 0.6;
    C_CRN *= 1.2;
	t = 0.0;
	udt = 0.0000025;
	end_time = 1.2;
	steps = end_time/udt;
	
	
	char c;
	double gx = 0.0, gy = 0.0, gz = 0.0;
	
	for (z = 0; z < Z; z++) {
		for (y = 0; y < Y; y++) {
			for (x = 0; x < X; x++) {     
				g[x][y][z] = 0;
				xx[x][y][z] = 1.0;
				yy[x][y][z] = 0.0;
				zz[x][y][z] = 0.0;
			}
		}
	}
	for(x = 2; x <= X-2; x++){
		for(y = 2; y <= Y-2; y++){
			for(z = 2; z <= Z-2; z++){
				g[x][y][z] = 1;
			}
		}
	}
	for(x = 2; x <= 11; x++){
		for(y = 2; y <= Y-2; y++){
			for(z = 2; z <= Z-2; z++){
				g[x][y][z] = 3;
			}
		}
	}
	


	

//*
	for (x = 1; x < X; x++) 
		for (y = 1; y < Y; y++) 
			for (z = 1; z < Z; z++)
				if (g[x][y][z] > 0)
					h[x][y][z] = 1;

	num = 1;      
	for (x = 1; x < X; x++) 
		for (y = 1; y < Y; y++) 
			for (z = 1; z < Z; z++)
				if (h[x][y][z] > 0){
					h[x][y][z] = num;
					num++;
				}

	for (x = 1; x < X; x++) 
		for (y = 1; y < Y; y++) 
			for (z = 1; z < Z; z++){
				gg[1] = h[x - 1][y - 1][z - 1];
				gg[2] = h[x - 1][y - 1][z];
				gg[3] = h[x - 1][y - 1][z + 1];
				gg[4] = h[x - 1][y][z - 1];
				gg[5] = h[x - 1][y][z];
				gg[6] = h[x - 1][y][z + 1];
				gg[7] = h[x - 1][y + 1][z - 1];
				gg[8] = h[x - 1][y + 1][z];
				gg[9] = h[x - 1][y + 1][z + 1];

				gg[10] = h[x][y - 1][z - 1];
				gg[11] = h[x][y - 1][z];
				gg[12] = h[x][y - 1][z + 1];
				gg[13] = h[x][y][z - 1];
				gg[14] = h[x][y][z + 1];
				gg[15] = h[x][y + 1][z - 1];
				gg[16] = h[x][y + 1][z];
				gg[17] = h[x][y + 1][z + 1];

				gg[18] = h[x + 1][y - 1][z - 1];
				gg[19] = h[x + 1][y - 1][z];
				gg[20] = h[x + 1][y - 1][z + 1];
				gg[21] = h[x + 1][y][z - 1];
				gg[22] = h[x + 1][y][z];
				gg[23] = h[x + 1][y][z + 1];
				gg[24] = h[x + 1][y + 1][z - 1];
				gg[25] = h[x + 1][y + 1][z];
				gg[26] = h[x + 1][y + 1][z + 1];

				isbound = 0;
				for(i = 1; i <= 26; i++){ 
					if (gg[i] > 0){
						gg[i] = 1;
						isbound++;
					} 
					else
						gg[i] = 0;
				}

				if (h[x][y][z] == 0 && isbound > 0) {
					ic = (gg[3]/root3) - (gg[1]/root3) + (gg[6]/root2) +
						(gg[9]/root3) - (gg[7]/root3) - (gg[4]/root2) +
						(gg[12]/root2) - (gg[10]/root2) + gg[14] +
						(gg[17]/root2) - (gg[15]/root2) - gg[13] +
						(gg[20]/root3) - (gg[18]/root3) + (gg[23]/root2) +
						(gg[26]/root3) - (gg[24]/root3) - (gg[21]/root2);

					ir = (gg[9]/root3) - (gg[2]/root2) - (gg[3]/root3) - 
						(gg[1]/root3) + (gg[8]/root2) + (gg[7]/root3) +
						(gg[17]/root2) - gg[11] - (gg[12]/root2) -
						(gg[10]/root2) + gg[16] + (gg[15]/root2) +
						(gg[26]/root3) - (gg[19]/root2) - (gg[20]/root3) -
						(gg[18]/root3) + (gg[25]/root2) + (gg[24]/root3);

					il = (gg[18]/root3) + (gg[19]/root2) + (gg[20]/root3) +
						(gg[21]/root2) + gg[22] + (gg[23]/root2) +
						(gg[24]/root3) + (gg[25]/root2) + (gg[26]/root3) -
						(gg[1]/root3) - (gg[2]/root2) - (gg[3]/root3) -
						(gg[4]/root2) - gg[5] - (gg[6]/root2) - 
						(gg[7]/root3) - (gg[8]/root2) - (gg[9]/root3);

					imax = fabs(ic);
					if (fabs(ir) > imax)
						imax = fabs(ir);
					if (fabs(il) > imax)
						imax = fabs(il);

					i = 0;
					j = 0; 
					k = 0;

					tflt = ir / fabs(imax);
					if (tflt <= 0.5 && tflt >= -0.5) 
						i = 0;
					else if (tflt > 0.5) 
						i = 1;
					else if (tflt < -0.5) 
						i = -1;

					tflt = ic / fabs(imax);
					if (tflt <= 0.5 && tflt >= -0.5) 
						j = 0;
					else if (tflt > 0.5) 
						j = 1;
					else if (tflt < -0.5) 
						j = -1;

					tflt = il / fabs(imax);
					if (tflt <= 0.5 && tflt >= -0.5)
						k = 0;
					else if (tflt > 0.5)
						k = 1;
					else if (tflt < -0.5)
						k = -1;

					if (imax == 0){
						i = 0;
						j = 0;
						k = 0;
					}
					if (h[x + k][y + i][z + j] > 0)    
						h[x][y][z] = -1 * h[x + k][y + i][z + j];   
					else
						h[x][y][z] = h[x + k][y + i][z + j];
				}
			}
//*/
	for (z = 0; z < Z; z++) {
		for (y = 0; y < Y; y++) {
			for (x = 0; x < X; x++) {
				if(g[x][y][z] == 0){ 
					Y0[x][y][z]  = 0.0;
					Y1[x][y][z]  = 0.0;
					Y2[x][y][z]  = 0.0;
					Y3[x][y][z]  = 0.0;
					Y10[x][y][z] = 0.0;
					Y11[x][y][z] = 0.0;
					Y14[x][y][z] = -81.2;
					Y15[x][y][z] = 0.0;
					Y16[x][y][z] = 0.0;
					Y17[x][y][z] = 0.0;
					Y18[x][y][z] = 0.0;
					Y22[x][y][z] = 0.0;
					Y24[x][y][z] = 0.0;
					Y25[x][y][z] = 0.0;
					Y26[x][y][z] = 0.0;
					Y27[x][y][z] = 0.0;
					Y28[x][y][z] = 0.0;
					Y29[x][y][z] = 0.0;
					Y30[x][y][z] = 0.0;
					Y31[x][y][z] = 0.0;
					Y32[x][y][z] = 0.0;
				} 
				else if ((g[x][y][z] == 1) || (g[x][y][z] == 3)){
					Y0[x][y][z] = 0.00; 
					Y1[x][y][z] = 1.00; 
					Y2[x][y][z] = 0.9992; 
					Y3[x][y][z] = 0.0; 
					Y10[x][y][z] = 1.488;  
					Y11[x][y][z] = 1.488;   
					Y12[x][y][z] = 139;    
					Y13[x][y][z] = 0.0001013; 
					Y14[x][y][z] = -81.2;   
					Y15[x][y][z] = 11.2;    
					Y16[x][y][z] = 0.000137; 
					Y17[x][y][z] = 0.775; 
					Y18[x][y][z] = 0.999837;
					Y22[x][y][z] = 0.0000329; 
					Y24[x][y][z] = 0.00254; 		
					Y25[x][y][z] = 0.01869; 
					Y26[x][y][z] = 0.999; 
					Y27[x][y][z] = 0.00496; 
					Y28[x][y][z] = 0.96500; 
					Y29[x][y][z] = 0.00291; 
					Y30[x][y][z] = 0.97800; 
					Y31[x][y][z] = 0.999; 
					Y32[x][y][z] = 0.0304; 	
				}
			}
		}
	}
	
	Cao = 1.8;   
	Ko = 5.4;  
	Nao = 140.0;   
	F = 96485.0;
	R2 = 8314.0;
	T = 310.0;   
	RTonF = R2*T/F;
	ACh_on = 1.0;   
	ACh = 0.0e-6;

	for (x = 1; x < X; x++)
    	for (y = 1; y < Y; y++)
          	for (z = 1; z < Z; z++)
      			if (g[x][y][z] > 0) {
					if((g[x][y][z] == 1)||(g[x][y][z] == 3)){
						D1 = D1A;
						D2 = D2A;
						DD = DDA;
					}
					dc[x][y][z][1] = D2 + (DD * xx[x][y][z] * xx[x][y][z]);
					dc[x][y][z][2] = DD * xx[x][y][z] * yy[x][y][z];
					dc[x][y][z][3] = DD * xx[x][y][z] * zz[x][y][z];
					dc[x][y][z][4] = DD * yy[x][y][z] * xx[x][y][z];
					dc[x][y][z][5] = D2 + (DD * yy[x][y][z] * yy[x][y][z]);
					dc[x][y][z][6] = DD * yy[x][y][z] * zz[x][y][z];
					dc[x][y][z][7] = DD * zz[x][y][z] * xx[x][y][z];
					dc[x][y][z][8] = DD * zz[x][y][z] * yy[x][y][z];
					dc[x][y][z][9] = D2 + (DD * zz[x][y][z] * zz[x][y][z]);
					df[x][y][z][1] = (xx[x + 1][y][z] - xx[x - 1][y][z]) / (2*dx);
					df[x][y][z][2] = (xx[x][y + 1][z] - xx[x][y - 1][z]) / (2*dx);
					df[x][y][z][3] = (xx[x][y][z + 1] - xx[x][y][z - 1]) / (2*dx);
					df[x][y][z][4] = (yy[x + 1][y][z] - yy[x - 1][y][z]) / (2*dx);
					df[x][y][z][5] = (yy[x][y + 1][z] - yy[x][y - 1][z]) / (2*dx);
					df[x][y][z][6] = (yy[x][y][z + 1] - yy[x][y][z - 1]) / (2*dx);
					df[x][y][z][7] = (zz[x + 1][y][z] - zz[x - 1][y][z]) / (2*dx);
					df[x][y][z][8] = (zz[x][y + 1][z] - zz[x][y - 1][z]) / (2*dx);
					df[x][y][z][9] = (zz[x][y][z + 1] - zz[x][y][z - 1]) / (2*dx); 
				}

	for (increment = 0; increment < steps; increment++){     
        if ((fabs(t - ptstim) < 0.00001)) {
			for (x = 1; x < X; x++){
    			for (y = 1; y < Y; y++){
          			for (z = 1; z < Z; z++){
      					if (g[x][y][z] == 3){
							Y14[x][y][z] = 20.0;
						}
					}
				}
			}
            ptstim = ptstim + counter;
        }
		for (x = 1; x < X; x++){
			for (y = 1; y < Y; y++){
				for (z = 1; z < Z; z++){ 
					if (h[x][y][z] < 0){
						for (i = -1; i <= 1; i++){
							for (j = -1; j <= 1; j++){
								for (k = -1; k <= 1; k++){
									if (h[x][y][z] == -h[x + i][y + j][z + k]){
										Y14[x][y][z] = Y14[x + i][y + j][z + k];
									}
								}
							}
						}
					}
				}
			}
		}
		
		for (x = 1; x < X; x++){
			for (y = 1; y < Y; y++){
				for (z = 1; z < Z; z++){
					if ((g[x][y][z] == 1)||(g[x][y][z] == 3)){ 
						i_rel = Y3[x][y][z];
						Ki =  Y12[x][y][z];
						V =  Y14[x][y][z];
						Nai = Y15[x][y][z];
						E_K = RTonF*log(Ko/Ki);
						if (V < -40) {
							alpha_h = 135.0*exp(-(80+V)/6.8);  
							beta_h = 3560.0*exp(0.079*V) + 310000000*exp(0.35*V);
						}
						else{
							alpha_h = 0.0;
							beta_h = 1/(0.00013*(1+exp(-(V+10.66)/11.1)));
						}
						h_infinity = alpha_h/(alpha_h+beta_h);
						tau_h = 1/(alpha_h+beta_h);
						hgate = Y28[x][y][z];
						dY28 = (h_infinity - hgate)/tau_h;
						if(V == -47.13)
							alpha_m = 3200.0;
						else
							alpha_m = 320.0*(V+47.13)/(1-exp(-0.1*(V+47.13)));
						beta_m = 80.0*exp(-V/11);
						m_infinity = alpha_m/(alpha_m+beta_m);
						tau_m = 1/(alpha_m+beta_m);
						m = Y29[x][y][z];
						dY29 = (m_infinity - m)/tau_m;
						if (V < -40) {
							alpha_j = (-127140000*exp(0.2444*V) - 0.03474*exp(-0.04391*V))*(V+37.78) / (1+exp(0.311*(V+79.23)));
							beta_j = 121.2*exp(-0.01052*V)/(1+exp(-0.1378*(V+40.14)));
						}
						else{
							alpha_j = 0.0;
							beta_j = 300.0*exp(-0.0000002535*V)/(1+exp(-0.1*(V+32)));
						}
						j_infinity = alpha_j/(alpha_j+beta_j);
						tau_j = 1/(alpha_j+beta_j);
						jgate = Y30[x][y][z];
						dY30 = (j_infinity - jgate)/tau_j;
						E_Na = RTonF*log(Nao/Nai);
						i_Na = gna*m*m*m*hgate*jgate*(V-E_Na);
						d = Y16[x][y][z];
						d_infinity = 1/(1+exp(-(V+10)/8));
						if(fabs(V+10)<1e-10)
							tau_d = 0.0045790/(1+exp(-(V+10)/6.24));
						else
							tau_d = (1-exp(-(V+10)/6.24)) / (35.0*(V+10)*(1+exp(-(V+10)/6.24)));
						dY16 = (d_infinity - d)/tau_d;
						fca = Y17[x][y][z];
						fca_infinity = 1/(1+Y13[x][y][z]/0.00035);
						tau_fca = 0.002;													
						dY17 = (fca_infinity - fca)/tau_fca;
						f = Y18[x][y][z];
						f_infinity = exp(-(V+28)/6.9)/(1+exp(-(V+28)/6.9));
						tau_f = 0.009/(0.0197*exp(-pow((0.0337*(V+10)),2))+0.02);
						dY18 = (f_infinity - f)/tau_f;
						i_CaL = 0.45*d*f*fca*gcalbar*(V-65);
						xr = Y22[x][y][z];
						if(fabs(V+14.1)<1e-10)
							alpha_xr = 1.5;
						else
							alpha_xr = 0.3*(V+14.1)/(1-(exp(-(V+14.1)/5)));       
						if(fabs(V-3.3328)<1e-10)
							beta_xr = 3.7836118e-1;
						else    
							beta_xr = 0.073898*(V - 3.3328)/(exp((V-3.3328)/5.1237) - 1);    
						tau_xr = 1/(alpha_xr+beta_xr);
						xr_infinity = 1/(1+exp(-(V+14.1)/6.5));
						dY22 = (xr_infinity - xr)/tau_xr;
						i_Kr = gkr*xr*(V-E_K)/(1+(exp((V+15)/22.4)));
						xs = Y25[x][y][z];
						if(fabs(V-19.9)<1e-10){
							alpha_xs = 0.68;
							beta_xs = 0.315;   
						}
						else{
							alpha_xs = 0.04*(V-19.9)/(1-exp(-(V-19.9)/17));
							beta_xs = 0.035*(V-19.9)/(exp((V-19.9)/9)-1);
						} 
						tau_xs = 0.5/(alpha_xs+beta_xs);
						xs_infinity = pow(1 + exp(-(V-19.9)/12.7),-0.5);                       
						dY25 = (xs_infinity - xs)/tau_xs;
						i_Ks = 2*gks*xs*xs*(V-E_K);
						i_K1 = 2*gk1*(V-E_K)/(1 + exp(0.07*(V+80)));
						gkur = gkurMultiplyer*(0.5 + 5.0/(1+exp(-(V-15)/13.0)));
						ui  = Y26[x][y][z];					
						alpha_ui = 1000/(21+exp(-(V-185)/28));
						beta_ui = 1000/(exp(-(V-158)/16));
						tau_ui = 1/(3*(alpha_ui+beta_ui));
						ui_infinity = 1/(1+exp((V-99.45)/27.48));
						dY26 = (ui_infinity - ui)/tau_ui;
						ua  = Y27[x][y][z];
						alpha_ua = 650.0/(exp(-(V+10)/8.5) + exp(-(V-30)/59.0));
						beta_ua = 650.0/(2.5 + exp((V+82)/17.0));
						tau_ua = 1/(3*(alpha_ua+beta_ua));
						ua_infinity = 1/(1+exp(-(V+30.3)/9.6));
						dY27 = (ua_infinity - ua)/tau_ua;
						i_Kur = 0.5*gkur*ua*ua*ua*ui*(V-E_K);
						igate = Y31[x][y][z];
						alpha_i = 1000/(18.53+exp((V+113.7)/10.95));
						beta_i = 1000/(35.56+exp(-(V+1.26)/7.44));                          
						tau_i = 1/(3*(alpha_i+beta_i));                          
						i_infinity = 1/(1+exp((V+43.1)/5.3));                          
						dY31 = (i_infinity - igate)/tau_i;
						a = Y32[x][y][z];
						alpha_a = 650.0/(exp(-(V+10)/8.5) + exp(-(V-30)/59));
						beta_a = 650.0/(2.5+exp((V+82)/17));                           
						tau_a = 1/(3*(alpha_a+beta_a));                           
						a_infinity = 1/(1+exp(-(V+20.47)/17.54));
						dY32 = (a_infinity - a)/tau_a;
						i_to = 0.35*gito*a*a*a*igate*(V-E_K);
						i_NaCa = 1.6*ibarpaca*(exp(0.35*V/RTonF)*pow(Nai,3)*Cao - exp((0.35-1)*V/RTonF)*pow(Nao,3)*Y13[x][y][z]) * ( 1/((pow(kmnancx,3)+pow(Nao,3))) * 1/(kmcancx+Cao) * 1/(1 + (0.1*exp((0.35-1)*V/RTonF))));
						i_NaK = ibarnak/(1 + 0.1245*exp(-0.1*V/RTonF) + 0.0365*(exp(Nao/67.3)-1)/7.0 * exp(-V/RTonF)) * 1/(1+pow((kmnai/Nai),1.5)) * Ko/(Ko+kmko);
						i_pCa = ibarpca*Y13[x][y][z]/(kmpca + Y13[x][y][z]);
						E_Ca = 0.5*RTonF*log(Cao/Y13[x][y][z]);
						i_bCa = gcab*(V-E_Ca);
						i_bK = gkb*(V-E_K);
						i_bNa = gnab*(V-E_Na);
						if(ACh > 0.0){
							yACh = Y24[x][y][z];
							alpha_yACh = 12.32/(1 + 0.0042/ACh) + 0.2475;
							beta_yACh = 10.0*exp(0.0133*(V + 40));
							tau_yACh = 1/(alpha_yACh + beta_yACh);
							yACh_infinity = alpha_yACh/(alpha_yACh + beta_yACh);
							dY24 = (yACh_infinity - yACh)/tau_yACh;
							i_KACh = ACh_on*gkach*yACh*(V - E_K)/(1+exp((V + 20)/20));
						}
						else{
							i_KACh = 0.0;
						}
						Na_ion_tot = i_Na + i_bNa + 3*i_NaK + 3*i_NaCa;
						K_ion_tot = i_Kr + i_Ks + i_K1 + i_bK - 2*i_NaK + i_to + i_Kur + i_KACh;
						Ca_ion_tot = i_CaL + i_bCa + i_pCa - 2*i_NaCa;
						i_tot = Na_ion_tot + K_ion_tot + Ca_ion_tot;
						dY14 = -i_tot/C_CRN;
						i_tr = (Y11[x][y][z] - Y10[x][y][z])/tautr;
						u = Y0[x][y][z];
						fn = 1e-6*vjsr*i_rel - 5e-8/F * (5*i_CaL - 2*i_NaCa);
						u_infinity = 1/(1 + exp(-(fn-3.4175e-13)/13.67e-16));
						dY0 = (u_infinity - u)/tau_u;
						v = Y1[x][y][z];
						tau_v = 0.00191 + 0.00209/(1 + exp(-(fn-3.4175e-13)/13.67e-16));
						v_infinity = 1 - 1/(1+exp(-(fn-6.835e-14)/13.67e-16));
						dY1 = (v_infinity - v)/tau_v;
						w = Y2[x][y][z];
						if(fabs(V-7.9)<1e-10)
							tau_w = 0.0602/1.3;
						else 
							tau_w = 0.006*(1-exp(-(V-7.9)/5)) / ((1+0.3*exp(-(V-7.9)/5))*(V-7.9));    
						w_infinity = 1 - 1/(1+exp(-(V-40)/17));
						dY2 = (w_infinity - w)/tau_w;
						i_rel = 30000*u*u*v*w*(Y10[x][y][z]-Y13[x][y][z]);
						dY10 = (i_tr-i_rel)/(1+(10*0.8/pow((Y10[x][y][z]+0.8),2)));
						dY15 = -Na_ion_tot/(vmyo*F);
						dY12 = -K_ion_tot/(vmyo*F);
						i_leak = 1.5*5*Y11[x][y][z]/15.0;
						i_up = 5*Y13[x][y][z]/(Y13[x][y][z] + 0.00092);
						dY11 = i_up - i_leak - i_tr*vjsr/vnsr;
						cmdn = cmdnbar*(Y13[x][y][z]/(Y13[x][y][z] + kmcmdn));
						b1Cai = -Ca_ion_tot/(2*F*vmyo) + ((i_leak-i_up)*vnsr + i_rel*vjsr)/vmyo;
						b2Cai = (1+(trpnbar*kmtrpn/pow((Y13[x][y][z]+kmtrpn),2))+(cmdnbar*kmcmdn/pow((Y13[x][y][z]+kmcmdn),2)));
						dY13 = b1Cai/b2Cai;	
					}
					if(g[x][y][z] > 0){
						dudx2  = (Y14[x - 1][y][z] + Y14[x + 1][y][z] - 2 * Y14[x][y][z]) / (dx*dx);          
						dudy2  = (Y14[x][y - 1][z] + Y14[x][y + 1][z] - 2 * Y14[x][y][z]) / (dx*dx);
						dudz2  = (Y14[x][y][z - 1] + Y14[x][y][z + 1] - 2 * Y14[x][y][z]) / (dx*dx);
						dudxdy = (Y14[x + 1][y + 1][z] + Y14[x - 1][y - 1][z] - Y14[x + 1][y - 1][z] - Y14[x - 1][y + 1][z])/(4*dx*dx);  
						dudxdz = (Y14[x + 1][y][z + 1] + Y14[x - 1][y][z - 1] - Y14[x + 1][y][z - 1] - Y14[x - 1][y][z + 1])/(4*dx*dx);
						dudydz = (Y14[x][y + 1][z + 1] + Y14[x][y - 1][z - 1] - Y14[x][y + 1][z - 1] - Y14[x][y - 1][z + 1])/(4*dx*dx);
						dudx   = (Y14[x + 1][y][z] - Y14[x - 1][y][z])/(2*dx);  
						dudy   = (Y14[x][y + 1][z] - Y14[x][y - 1][z])/(2*dx);
						dudz   = (Y14[x][y][z + 1] - Y14[x][y][z - 1])/(2*dx);
						D1 = D1A;
						D2 = D2A;
						DD = DDA;
						du= dc[x][y][z][1]*dudx2  + dudx*DD*(xx[x][y][z]*df[x][y][z][1] + xx[x][y][z]*df[x][y][z][1]) +
							dc[x][y][z][2]*dudxdy + dudy*DD*(xx[x][y][z]*df[x][y][z][4] + yy[x][y][z]*df[x][y][z][1]) +
							dc[x][y][z][3]*dudxdz + dudz*DD*(xx[x][y][z]*df[x][y][z][7] + zz[x][y][z]*df[x][y][z][1]) +
							dc[x][y][z][4]*dudxdy + dudx*DD*(yy[x][y][z]*df[x][y][z][2] + xx[x][y][z]*df[x][y][z][5]) +
							dc[x][y][z][5]*dudy2  + dudy*DD*(yy[x][y][z]*df[x][y][z][5] + yy[x][y][z]*df[x][y][z][5]) +
							dc[x][y][z][6]*dudydz + dudz*DD*(yy[x][y][z]*df[x][y][z][8] + zz[x][y][z]*df[x][y][z][5]) +
							dc[x][y][z][7]*dudxdz + dudx*DD*(zz[x][y][z]*df[x][y][z][3] + xx[x][y][z]*df[x][y][z][9]) +
							dc[x][y][z][8]*dudydz + dudy*DD*(zz[x][y][z]*df[x][y][z][6] + yy[x][y][z]*df[x][y][z][9]) +
							dc[x][y][z][9]*dudz2  + dudz*DD*(zz[x][y][z]*df[x][y][z][9] + zz[x][y][z]*df[x][y][z][9]);
						Y0[x][y][z] = Y0[x][y][z] + udt*dY0;
						Y1[x][y][z] = Y1[x][y][z] + udt*dY1;
						Y2[x][y][z] = Y2[x][y][z] + udt*dY2;
						Y3[x][y][z] = i_rel;					
						Y10[x][y][z] = Y10[x][y][z] + udt*dY10;
						Y11[x][y][z] = Y11[x][y][z] + udt*dY11;
						Y12[x][y][z] = Y12[x][y][z] + udt*dY12;
						Y13[x][y][z] = Y13[x][y][z] + udt*dY13;
						Y14[x][y][z] = Y14[x][y][z] + udt*(du/1.2 + dY14);
						Y15[x][y][z] = Y15[x][y][z] + udt*dY15;
						Y16[x][y][z] = Y16[x][y][z] + udt*dY16;
						Y17[x][y][z] = Y17[x][y][z] + udt*dY17;
						Y18[x][y][z] = Y18[x][y][z] + udt*dY18;						
						Y22[x][y][z] = Y22[x][y][z] + udt*dY22;
						Y24[x][y][z] = Y24[x][y][z] + udt*dY24;
						Y25[x][y][z] = Y25[x][y][z] + udt*dY25;
						Y26[x][y][z] = Y26[x][y][z] + udt*dY26;
						Y27[x][y][z] = Y27[x][y][z] + udt*dY27;
						Y28[x][y][z] = Y28[x][y][z] + udt*dY28;
						Y29[x][y][z] = Y29[x][y][z] + udt*dY29;
						Y30[x][y][z] = Y30[x][y][z] + udt*dY30;
						Y31[x][y][z] = Y31[x][y][z] + udt*dY31;
						Y32[x][y][z] = Y32[x][y][z] + udt*dY32;
					}
				}
			}
		}
		t = t+udt;
		if ((increment % incrementDivisor == 0)&&(t>0.0)) {
			str = malloc(400 * sizeof(char));
			sprintf(str, "test/aa%.5ld.vtk", increment / incrementDivisor);
			printf("File number %ld, time = %f s, %.2f%% complete\n", increment / incrementDivisor, t, 100*t/end_time);
			out = fopen(str, "wt");
			for (z = 0; z < Z; z++) {
				for (y = 0; y < Y; y++) {
					for (x = 0; x < X; x++) {
						fprintf(out, "%2.2f ", Y14[x][y][z]);
					}
					fprintf(out, "\n");
				}
				fprintf(out, "\n");
			}
			fclose(out);
			free(str);
		}
	}
	return 0;
}