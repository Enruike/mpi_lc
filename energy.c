#include "finite.h"

void free_energy(){

	//Cambiaremos la dimensión de la variable a 3 por ser 2 regiones y una total.
	//double p_en_ldg = 0;
	double p_en_ldg[3]  = { 0. };
	double p_en_surf[2] = { 0. };
    double p_en_el[5] = { 0. };

	//Region interna
	double p_en_el_in[5] = { 0. };
	//Region externa
	double p_en_el_out[5] = { 0. };

	energy_ldg(p_en_ldg);                                                                      
	energy_el(p_en_el, p_en_el_in, p_en_el_out);                                                                           
	energy_surf(p_en_surf);                                                                    
	MPI_Barrier(MPI_COMM_WORLD);                                                                  
	MPI_Win_fence(0, win);                                                                        

	//Sum up the energy fraction in different processors and calculate dE
	MPI_Allreduce(&p_en_ldg, &en_ldg, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);                 
	MPI_Allreduce(&p_en_el, &en_el, 5, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);                     
	MPI_Allreduce(&p_en_surf, &en_surf, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);   
	//Agregamos más variables que sumarán los valores de las regiones interna y externa.
	MPI_Allreduce(&p_en_el_in, &en_el_in, 5, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(&p_en_el_out, &en_el_out, 5, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	
	en_tot = en_ldg[0] + en_el[0] + en_el[1] + en_el[2] + en_el[3] + en_el[4] + en_surf[0] + en_surf[1];                                  
	dE = en_tot - el_old;                                                             
    el_old = en_tot;
	
	if(DoubleU){
		if(myid == root && cycle % check_every == 0){
			printf("En_LDG: %lf, En_L1: %lf, ", en_ldg[0], en_el[0]);
			if(L2 != 0){
				printf("L2: %lf, ", en_el[1]);
			}
			if(chiral != 0) printf("En_Chiral: %lf ", en_el[4]);
			printf("En_Surf1: %lf, Par_Surf: %lf, Cycle: %d\n", en_surf[0], en_surf[1], cycle);
			if(geo != 10){
				printf("\t\t ~Inner Region Energy~\n");
			}
			else{
				printf("\t\t ~Energy~\n");
			}
			printf("LdG: %lf, L1: %lf, ", en_ldg[1], en_el_in[0]);
			if(L2 != 0){
				printf("L2: %lf, ", en_el_in[1]);
			}
			printf("Chiral: %lf\n", en_el_in[4]);
			if(geo != 10){
				printf("\t\t ~Outer Region Energy~\n");
			}
			else{
				printf("\t ~Energy Around the Interface~\n");
			}
			printf("LdG: %lf, L1: %lf, ", en_ldg[2], en_el_out[0]);
			if(L2 != 0){
				printf("L2: %lf, ", en_el_out[1]);
			}
			printf("Chiral: %lf\n", en_el_out[4]);
			printf("Total Energy: %lf \n", en_tot);		
			printf("dE: %lf \n\n", dE);		
		}
	}
	else{
		if(myid == root && cycle % check_every == 0){
			printf("En_LDG: %lf, En_L1: %lf, En_L2: %lf, En_L3: %lf, En_L4: %lf, En_Chiral: %lf, En_Surf1: %lf, En_Surf2: %lf Cycle: %d\n", \
			en_ldg[0], en_el[0], en_el[1], en_el[2], en_el[3], en_el[4], en_surf[0], en_surf[1], cycle);
			printf("Total Energy: %lf \n", en_tot);		
			printf("dE: %lf \n\n", dE);		
		}
	}

}	

//Landau de-Gennes energy for bulk points.
void energy_ldg(double* ans){
	double trace2 = 0.;
	double trace3 = 0.;
	double Qin[6] = { 0. };

	if(DoubleU && geo == -44 && flag_2 == false){
		for(int i = 0; i < length; i++)	{
			if(sign[i] == 0 || sign[i] == 1){
				
				Qin[0] = q[i * 6 + 0];
				Qin[1] = q[i * 6 + 1];
				Qin[2] = q[i * 6 + 2];
				Qin[3] = q[i * 6 + 3];
				Qin[4] = q[i * 6 + 4];
				Qin[5] = q[i * 6 + 5];

				trace2 = trqq(Qin);
				trace3 = trqqq(Qin);
				
				ans[0] += 0.5 * (1. - U2 / 3.) * trace2 - U2 / 3. * trace3 + U2 * 0.25 * trace2 * trace2;
				
			}
		}
	}
	else if(DoubleU){
		for (int i = 0; i < length; i ++){				

			if(sign[i] == 0 || sign[i] == 1){
				
				Qin[0] = q[i * 6 + 0];
				Qin[1] = q[i * 6 + 1];
				Qin[2] = q[i * 6 + 2];
				Qin[3] = q[i * 6 + 3];
				Qin[4] = q[i * 6 + 4];
				Qin[5] = q[i * 6 + 5];

				trace2 = trqq(Qin);
				trace3 = trqqq(Qin);

				if(bulktype_MPI[i] == 1){
					ans[0] += 0.5 * (1. - U / 3.) * trace2 - U / 3. * trace3 + U * 0.25 * trace2 * trace2;
					ans[1] += 0.5 * (1. - U / 3.) * trace2 - U / 3. * trace3 + U * 0.25 * trace2 * trace2;
				}
				
				else if(bulktype_MPI[i] == 2 || bulktype_MPI[i] == 3){
					ans[0] += 0.5 * (1. - U2 / 3.) * trace2 - U2 / 3. * trace3 + U2 * 0.25 * trace2 * trace2;
					ans[2] += 0.5 * (1. - U2 / 3.) * trace2 - U2 / 3. * trace3 + U2 * 0.25 * trace2 * trace2;
				}
			}
		}
	}
	else{

		for (int i = 0; i < length; i ++){				

			if(sign[i] == 0 || sign[i] == 1){

				Qin[0] = q[i * 6 + 0];
				Qin[1] = q[i * 6 + 1];
				Qin[2] = q[i * 6 + 2];
				Qin[3] = q[i * 6 + 3];
				Qin[4] = q[i * 6 + 4];
				Qin[5] = q[i * 6 + 5];	
			
				trace2 = trqq(Qin);
				trace3 = trqqq(Qin);

				ans[0] += 0.5 * (1. - U / 3.) * trace2 - U / 3. * trace3 + U * 0.25 * trace2 * trace2;
				//if (i % 10 == 0 && cycle % 50 == 0 && sign[i] == 1) printf("traceqq = %lf  trqqq = %lf i = %d sign[%d] = %d \n", trqq(Qin), trqqq(Qin), i, i, sign[i]);
			}
		}
	}
	
	if(DoubleU){
		ans[0] *= dV;
		ans[1] *= dVi;
		ans[2] *= dVo;
	}
	else{
		ans[0] *= dV;
	}
}

//Elastic energy for bulk points
void energy_el(double* ans, double* ans_in, double* ans_out){
	int i, n, j, k, l;
	double dQ[3][6];
	double Qin[6] = { 0. };
	double vec[3] = { 0. };
	int ref = length * myid;
	int xm, xp, ym, yp, zm, zp;
	for (i = 0; i < 3; i ++){
		for(j = 0; j < 6; j ++){
			dQ[i][j] = 0;
		}
	}
	//for(i = 0; i < 5; i ++)	ans[i] = 0;

	for (i = 0; i < length; i ++){
		if(sign[i] == 0 || sign[i] == 1){
			for(n = 0; n < 6; n ++)	Qin[n] = q[i * 6 + n];
			xm = neigb[i * 6 + 0] - ref;
			xp = neigb[i * 6 + 1] - ref;
			ym = neigb[i * 6 + 2] - ref;
			yp = neigb[i * 6 + 3] - ref;
			zm = neigb[i * 6 + 4] - ref;
			zp = neigb[i * 6 + 5] - ref;
			for (n = 0; n < 6; n ++) {
				//dQ is the first derivative with second order approximation
				//first index for direction: 0-x; 1-y; 2-z;
				//second index for qtensor index;
				dQ[0][n] = (q[xp * 6 + n] - q[xm * 6 + n]) * 0.5 * idx;
				dQ[1][n] = (q[yp * 6 + n] - q[ym * 6 + n]) * 0.5 * idy;
				dQ[2][n] = (q[zp * 6 + n] - q[zm * 6 + n]) * 0.5 * idz;
			}
			ans[0] += trqq(dQ[0])+trqq(dQ[1])+trqq(dQ[2]);
			if(DoubleU){
				if(bulktype_MPI[i] == 1){
					ans_in[0] += trqq(dQ[0])+trqq(dQ[1])+trqq(dQ[2]);
				}
				else if(bulktype_MPI[i] == 2 || bulktype_MPI[i] == 3){
					ans_out[0] += trqq(dQ[0])+trqq(dQ[1])+trqq(dQ[2]);
				}
			}
			if (L2 != 0){
				vec[0] = dQ[0][0];
				vec[1] = dQ[1][1];
				vec[2] = dQ[2][2];
				ans[1] += matr_mult(vec);

				if(DoubleU){
					if(bulktype_MPI[i] == 1){
						ans_in[1] += matr_mult(vec);
					}
					else if(bulktype_MPI[i] == 2 || bulktype_MPI[i] == 3){
						ans_out[1] += matr_mult(vec);
					}
				}

				vec[0] = dQ[0][1];
				vec[1] = dQ[1][3];
				vec[2] = dQ[2][4];
				ans[1] += matr_mult(vec);

				if(DoubleU){
					if(bulktype_MPI[i] == 1){
						ans_in[1] += matr_mult(vec);
					}
					else if(bulktype_MPI[i] == 2 || bulktype_MPI[i] == 3){
						ans_out[1] += matr_mult(vec);
					}
				}

				vec[0] = dQ[0][2];
				vec[1] = dQ[1][4];
				vec[2] = dQ[2][5];
				ans[1] += matr_mult(vec);

				if(DoubleU){
					if(bulktype_MPI[i] == 1){
						ans_in[1] += matr_mult(vec);
					}
					else if(bulktype_MPI[i] == 2 || bulktype_MPI[i] == 3){
						ans_out[1] += matr_mult(vec);
					}
				}
				
			}
			if (L3 != 0){
				ans[2] += Qin[0] * trqq(dQ[0]) + Qin[3] * trqq(dQ[1]) + Qin[5] * trqq(dQ[2]) + 2 * Qin[1] * q_mult(dQ[0], dQ[1]) + 2 * Qin[2] * q_mult(dQ[0], dQ[2]) + 2 * Qin[4] * q_mult(dQ[1], dQ[2]); 
				if(DoubleU){
					if(bulktype_MPI[i] == 1){
						ans_in[2] += Qin[0] * trqq(dQ[0]) + Qin[3] * trqq(dQ[1]) + Qin[5] * trqq(dQ[2]) + 2 * Qin[1] * q_mult(dQ[0], dQ[1]) + 2 * Qin[2] * q_mult(dQ[0], dQ[2]) + 2 * Qin[4] * q_mult(dQ[1], dQ[2]);
					}
					else if(bulktype_MPI[i] == 2 || bulktype_MPI[i] == 3){
						ans_out[2] += Qin[0] * trqq(dQ[0]) + Qin[3] * trqq(dQ[1]) + Qin[5] * trqq(dQ[2]) + 2 * Qin[1] * q_mult(dQ[0], dQ[1]) + 2 * Qin[2] * q_mult(dQ[0], dQ[2]) + 2 * Qin[4] * q_mult(dQ[1], dQ[2]);
					}
				}
			
			}
			if (L4 != 0){
				ans[3] += dQ[0][0] * dQ[0][0] + dQ[1][1] * dQ[1][1] + dQ[2][2] * dQ[2][2] + dQ[0][1] * dQ[0][1] + dQ[1][3] * dQ[1][3] + dQ[2][4] * dQ[2][4] + dQ[0][2] * dQ[0][2] + dQ[1][4] * dQ[1][4] + dQ[2][5] * dQ[2][5] + 2 * (dQ[0][1] * dQ[1][0] + dQ[0][2] * dQ[2][0] + dQ[1][2] * dQ[2][1] + dQ[0][3] * dQ[1][1] + dQ[0][4] * dQ[2][1] + dQ[1][4] * dQ[2][3] 	+ dQ[1][2] * dQ[0][4] + dQ[0][5] * dQ[2][2] + dQ[1][5] * dQ[2][4]);
				if(DoubleU){
					if(bulktype_MPI[i] == 1){
						ans_in[3] += dQ[0][0] * dQ[0][0] + dQ[1][1] * dQ[1][1] + dQ[2][2] * dQ[2][2] + dQ[0][1] * dQ[0][1] + dQ[1][3] * dQ[1][3] + dQ[2][4] * dQ[2][4] + dQ[0][2] * dQ[0][2] + dQ[1][4] * dQ[1][4] + dQ[2][5] * dQ[2][5] + 2 * (dQ[0][1] * dQ[1][0] + dQ[0][2] * dQ[2][0] + dQ[1][2] * dQ[2][1] + dQ[0][3] * dQ[1][1] + dQ[0][4] * dQ[2][1] + dQ[1][4] * dQ[2][3] 	+ dQ[1][2] * dQ[0][4] + dQ[0][5] * dQ[2][2] + dQ[1][5] * dQ[2][4]);
					}
					else if(bulktype_MPI[i] == 2 || bulktype_MPI[i] == 3){
						ans_out[3] += dQ[0][0] * dQ[0][0] + dQ[1][1] * dQ[1][1] + dQ[2][2] * dQ[2][2] + dQ[0][1] * dQ[0][1] + dQ[1][3] * dQ[1][3] + dQ[2][4] * dQ[2][4] + dQ[0][2] * dQ[0][2] + dQ[1][4] * dQ[1][4] + dQ[2][5] * dQ[2][5] + 2 * (dQ[0][1] * dQ[1][0] + dQ[0][2] * dQ[2][0] + dQ[1][2] * dQ[2][1] + dQ[0][3] * dQ[1][1] + dQ[0][4] * dQ[2][1] + dQ[1][4] * dQ[2][3] 	+ dQ[1][2] * dQ[0][4] + dQ[0][5] * dQ[2][2] + dQ[1][5] * dQ[2][4]);
					}
				}
			
			}
			if(chiral == 1){
				//Chiral elastic energy
				ans[4] += Qin[0] * dQ[1][2] + Qin[1] * dQ[1][4] + Qin[2] * dQ[1][5] + Qin[1] * dQ[2][0] + Qin[3] * dQ[2][1] + Qin[4] * dQ[2][2]  + Qin[2] * dQ[0][1] + Qin[4] * dQ[0][3] + Qin[5] * dQ[0][4]  - Qin[0] * dQ[2][1] - Qin[1] * dQ[2][3] - Qin[2] * dQ[2][4]  - Qin[2] * dQ[1][0] - Qin[4] * dQ[1][1] - Qin[5] * dQ[1][2] - Qin[1] * dQ[0][2] - Qin[3] * dQ[0][4] - Qin[4] * dQ[0][5];
				
				if(DoubleU){
					if(bulktype_MPI[i] == 1){
						ans_in[4] += Qin[0] * dQ[1][2] + Qin[1] * dQ[1][4] + Qin[2] * dQ[1][5] + Qin[1] * dQ[2][0] + Qin[3] * dQ[2][1] + Qin[4] * dQ[2][2]  + Qin[2] * dQ[0][1] + Qin[4] * dQ[0][3] + Qin[5] * dQ[0][4]  - Qin[0] * dQ[2][1] - Qin[1] * dQ[2][3] - Qin[2] * dQ[2][4]  - Qin[2] * dQ[1][0] - Qin[4] * dQ[1][1] - Qin[5] * dQ[1][2] - Qin[1] * dQ[0][2] - Qin[3] * dQ[0][4] - Qin[4] * dQ[0][5];					}
					else if(bulktype_MPI[i] == 2 || bulktype_MPI[i] == 3){
						ans_out[4] += Qin[0] * dQ[1][2] + Qin[1] * dQ[1][4] + Qin[2] * dQ[1][5] + Qin[1] * dQ[2][0] + Qin[3] * dQ[2][1] + Qin[4] * dQ[2][2]  + Qin[2] * dQ[0][1] + Qin[4] * dQ[0][3] + Qin[5] * dQ[0][4]  - Qin[0] * dQ[2][1] - Qin[1] * dQ[2][3] - Qin[2] * dQ[2][4]  - Qin[2] * dQ[1][0] - Qin[4] * dQ[1][1] - Qin[5] * dQ[1][2] - Qin[1] * dQ[0][2] - Qin[3] * dQ[0][4] - Qin[4] * dQ[0][5];					}
				}

			}
		}
	}
	if(DoubleU){
		ans_in[0] *= 0.5 * dVi * L1;
		ans_in[1] *= 0.5 * dVi * L2;
		ans_in[2] *= 0.5 * dVi * L3;
		ans_in[3] *= 0.5 * dVi * L4;
		ans_in[4] *= dVi * (double)chiral * 2. * L1 * qch;

		ans_out[0] *= 0.5 * dVo * L1;
		ans_out[1] *= 0.5 * dVo * L2;
		ans_out[2] *= 0.5 * dVo * L3;
		ans_out[3] *= 0.5 * dVo * L4;
		ans_out[4] *= dVo * (double)chiral * 2. * L1 * qch;
	}
	ans[0] *= 0.5 * dV * L1;
	ans[1] *= 0.5 * dV * L2;
	ans[2] *= 0.5 * dV * L3;
	ans[3] *= 0.5 * dV * L4;
	ans[4] *= dV * (double)chiral * 2. * L1 * qch;
}

void energy_surf(double* ans){
	int i, n, nb;
	double Qdiff[6] = {0};
	double Qin[6] = {0};
	double loc_nu[3] = {0};
	int degen = 0, inf = 1;
	double Wstr = 0;
	double dA = 0;
	bool npboundary = true;
	nb = 0;
	for (i = 0; i < length; i++){
		if(sign[i] >= 2 && sign[i] <= 8 || sign[i] == 12 || sign[i] == 13 || (sign[i] >= 20 && sign[i] <= 23)){
			//for channel boundary
			if(sign[i] == 2 || sign[i] == 3){
				degen = degenerate;
				inf = infinite;
				Wstr = W;
				dA = dAdrop;
				npboundary = false;
			}
			//Para superficie inferior degenerada.
			else if(sign[i] == 12 || sign[i] == 13){
				degen = 1;
				inf = infinite;
				Wstr = W;
				dA = dAdrop;
				npboundary = false;
			}

			//for nanoparticle boundary
			else if(sign[i] == 4 || sign[i] == 5){
				degen = 0;
				inf = 0;
				Wstr = Wp;
				dA = dApart;
				npboundary = true;
			}
			else if(sign[i] == 6 || sign[i] == 7){
				degen = 1;
				inf = 0;
				Wstr = Wp;
				dA = dApart;
				npboundary = true;
			}
			else if(sign[i] >= 20 && sign[i] <= 23){
				degen = 2;
				inf = 0;
				Wstr = Wp;
				dA = dApart;
				npboundary = true;
			}
			else if(sign[i] == 8){
				continue;
			}
			else{
				printf("Error in energy_surf.\n");
			}

			if(Wstr != 0 && inf != 1){
				//printf("Test. sign[i] = %d\n", sign[i]);
				if(degen != 0){
					for(n = 0; n < 6; n ++)	Qin[n] = q[i * 6 + n];
					for(n = 0; n < 3; n ++)	loc_nu[n] = nu_p[nb * 3 + n];
					
					if(degen == 1) en_degen(Qin, loc_nu, Qdiff);
					else if(degen == 2) en_conic(Qin, loc_nu, Qdiff);

					if(npboundary){
						ans[1] += Wstr * trqq(Qdiff) * dApart;
					}	
					else{
						ans[0] += Wstr * trqq(Qdiff) * dAdrop;
					}
					
				}
				
				else if(degen == 0 && inf == 0){
					for(n = 0; n < 6; n ++){
						Qdiff[n] = q[i * 6 + n] - qo_p[nb * 6 + n];
					}
					if(npboundary){
						ans[1] += Wstr * trqq(Qdiff) * dApart;
					}	
					else{
						ans[0] += Wstr * trqq(Qdiff) * dAdrop;
					}
				}
			}
			nb ++;
		}
	}
}

void en_degen(double* Qin, double* loc_nu, double* Qdiff){
	int i, n, j, l, m;
	double Qtemp[3][3];
	double ptemp[3][3];
	double Qp[3][3];
	double third = 1.0 / 3;
	Qtemp[0][0] = Qin[0] + third * S;
	Qtemp[0][1] = Qtemp[1][0] = Qin[1];
	Qtemp[0][2] = Qtemp[2][0] = Qin[2];
	Qtemp[1][1] = Qin[3] + third * S;
	Qtemp[1][2] = Qtemp[2][1] = Qin[4];
	Qtemp[2][2] = Qin[5] + third * S;
	for(i = 0; i < 3; i++){
		for(j = 0; j < 3; j++){
			if(i == j) ptemp[i][j] = 1 - loc_nu[i] * loc_nu[j];
			else ptemp[i][j] = - loc_nu[i] * loc_nu[j];
		}
	}
	for(i = 0; i < 3; i++){
		for(j = 0; j < 3; j++){
			Qp[i][j] = 0;
			for(l = 0; l < 3; l++){
				for(m = 0; m < 3; m++){
					Qp[i][j] += ptemp[i][l]*Qtemp[l][m]*ptemp[m][j];
				}
			}
		}
	}
	Qdiff[0] = Qtemp[0][0] - Qp[0][0];
	Qdiff[1] = Qtemp[0][1] - Qp[0][1];
	Qdiff[2] = Qtemp[0][2] - Qp[0][2];
	Qdiff[3] = Qtemp[1][1] - Qp[1][1];
	Qdiff[4] = Qtemp[1][2] - Qp[1][2];
	Qdiff[5] = Qtemp[2][2] - Qp[2][2];
}

void en_conic(double* Qin, double* loc_nu, double* Qdiff){
	double Qtemp[3][3];
	double ptemp[3][3];
	double Qp[3][3];
	double third = 1. / 3.;
	double cosTiltAngle;
	double cosTiltAngleSq;
	
	Qtemp[0][0] = Qin[0] + third * S;
	Qtemp[0][1] = Qtemp[1][0] = Qin[1];
	Qtemp[0][2] = Qtemp[2][0] = Qin[2];
	Qtemp[1][1] = Qin[3] + third * S;
	Qtemp[1][2] = Qtemp[2][1] = Qin[4];
	Qtemp[2][2] = Qin[5] + third * S;

	for(int i = 0; i < 3; i++){
		for(int j = 0; j < 3; j++){
			ptemp[i][j] = loc_nu[i] * loc_nu[j];
		}
	}

	for(int i = 0; i < 3; i++){
		for(int j = 0; j < 3; j++){
			Qp[i][j] = 0;
			for(int l = 0; l < 3; l++){
				for(int m = 0; m < 3; m++){
					Qp[i][j] += ptemp[i][l] * Qtemp[l][m] * ptemp[m][j];
				}
			}
		}
	}
	
	cosTiltAngle = cos(tiltAngle / 180.0 * M_PI);
	cosTiltAngleSq = pow(cosTiltAngle, 2);
	
	Qdiff[0] =  Qp[0][0] - cosTiltAngleSq * S * ptemp[0][0];
	Qdiff[1] =  Qp[0][1] - cosTiltAngleSq * S * ptemp[0][1];
	Qdiff[2] =  Qp[0][2] - cosTiltAngleSq * S * ptemp[0][2];
	Qdiff[3] =  Qp[1][1] - cosTiltAngleSq * S * ptemp[1][1];
	Qdiff[4] =  Qp[1][2] - cosTiltAngleSq * S * ptemp[1][2];
	Qdiff[5] =  Qp[2][2] - cosTiltAngleSq * S * ptemp[2][2];
}