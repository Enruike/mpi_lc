#include "finite.h"

bool initial_ellip(){
	int l;
	int nb, nd, np;
	int i, j, k, n, m;
	double dis, x, y, z;
	double x1, y1, z1;

	//Center of the grid
	rx = lrint(Nx / 2);
	ry = lrint(Ny / 2);
	rz = lrint(Nz / 2);

	int xm, xp, ym, yp, zm, zp;
	int nsurf;
	bool *ndrop;
	int *indx;
	double **pos;
	int count1;
	
	int bulkin = 0;
	int bulkout = 0;

	dx = Lx / (double)(Nx - 1);
	dy = Ly / (double)(Ny - 1);
	dz = Lz / (double)(Nz - 1);

	double Rx = Lx / 2. - 2.;
	double Ry = Ly / 2. - 2.;
	double Rz = Lz / 2. - 2.;

	idx = 1. / dx;
	idy = 1. / dy;
	idz = 1. / dz;

	iddx = idx * idx;
	iddy = idy * idy;
	iddz = idz * idz;

	bulk = 0;
	surf = 0;
	nsurf = 0;
	tot = Nx * Ny * Nz;

	//allocate drop and boundary
	ndrop = (bool*)malloc(tot * sizeof(bool));
	nboundary = (bool*)malloc(tot * sizeof(bool));
	drop = (bool*)malloc(tot * sizeof(bool));
	boundary = (bool*)malloc(tot * sizeof(bool));
	indx = (int*)malloc(tot * sizeof(int));
	init_bulktype = (int*)malloc(tot * sizeof(int));

	//Crearemos dos índices para diferenciar el bulk del otro con una U distinta.

	for(int l = 0; l < tot; l ++){
		drop[l] = false;
		boundary[l] = false;
		ndrop[l] = false;
		nboundary[l] = false;
		indx[l] = -1;

		//Inicializamos nuestros nuevos índices

		init_bulktype[l] = 0;

	}

	l = 0;

	//Agregué un if que determina si la condición de doble U está activa. 
	if(DoubleU){

		//define the droplet 
		//Definiremos la gota oblatada completa para después generar la gota interna.
		for(int k = 0; k < Nz; k++){
			for (int j = 0; j < Ny; j++){
				for (int i = 0; i < Nx; i++){

					x = (double)(i - rx) * dx;
					y = (double)(j - ry) * dy;
					z = (double)(k - rz) * dz;

					if (( (x * x) / ((Rx + 0.5) * (Rx + 0.5)) + (y * y) / ((Ry + 0.5) * (Ry + 0.5)) + (z * z) / ((Rz + 0.5) * (Rz + 0.5)) ) <= 1 ){
						drop[l] = true;
						bulk++;

						//Añadimos un if anidado para determinar el bulk externo.
						if(( (x * x) / ((iRx + 0.5) * (iRx + 0.5)) + (y * y) / ((iRy + 0.5) * (iRy + 0.5)) + (z * z) / ((iRz + 0.5) * (iRz + 0.5)) ) <= 1){
							init_bulktype[l] = 1; //bulk interno.
							bulkin++;
						}
						else{
							init_bulktype[l] = 2; //para bulk externo.
							bulkout++;
						}

					}
					l ++;
				}
			}
		}
	}
	else{
		//define the droplet 
		//Definiremos la gota oblatada completa para después generar la gota interna.
		for(int k = 0; k < Nz; k++){
			for (int j = 0; j < Ny; j++){
				for (int i = 0; i < Nx; i++){
					x = (double)(i - rx) * dx;
					y = (double)(j - ry) * dy;
					z = (double)(k - rz) * dz;
					if(( (x * x) / ((Rx + 0.5) * (Rx + 0.5)) + (y * y) / ((Ry + 0.5) * (Ry + 0.5)) + (z * z) / ((Rz + 0.5) * (Rz + 0.5)) ) <= 1){
						drop[l] = true;
						bulk ++;
					}
					l ++;
				}
			}
		}
	}

	droplet = bulk;

	//define boundary
	l = 0;
	for(k = 0; k < Nz; k++){
		for (j = 0; j < Ny; j++){
			for (i = 0; i < Nx; i++){
				if(drop[l]){
					xm = i - 1 + j * Nx + k * Nx * Ny;
					xp = i + 1 + j * Nx + k * Nx * Ny;
					ym = i + (j - 1) * Nx + k * Nx * Ny;
					yp = i + (j + 1) * Nx + k * Nx * Ny;
					zm = i + j * Nx + (k - 1) * Nx * Ny;
					zp = i + j * Nx + (k + 1) * Nx * Ny;
					if(!drop[xm] || !drop[xp] || !drop[ym] || !drop[yp] || !drop[zm] || !drop[zp]){
						boundary[l] = true;
						surf ++;
					}
				}
				l ++;
			}
		}
	}
	bulk -= surf;
	for(int l = 0; l < tot; l ++){
		if(boundary[l])		drop[l] = false;
	}


	//Initialize particle if Np > 0
	if(Np != 0){
		pos = (double**)malloc(Np * sizeof(double*));
		for(i = 0; i < Np; i ++){
			pos[i] = (double*)malloc(4 * sizeof(double));
			for(n = 0; n < 4; n ++){
				pos[i][n] = 0;
			}
		}
		if(!read_nppos(pos)) return false;

		l = 0;
		for(k = 0; k < Nz; k++){
			for (j = 0; j < Ny; j++){
				for (i = 0; i < Nx; i++){
					if(boundary[l] || drop[l]){
						for(m = 0; m < Np; m ++){
							x = i - pos[m][0];
							y = j - pos[m][1];
							z = k - pos[m][2];
							dis = x*x+y*y+z*z;
							if (dis < (Rp - 0.5) * (Rp - 0.5)){
								ndrop[l] = true;
								if(drop[l]){
									bulk --;
									drop[l] = false;
								}
								else{	
									boundary[l] = false;
									surf --;
								}
							}
						}						
					}
					l ++;
				}
			}
		}
		l = 0;
		for(k = 0; k < Nz; k++){
			for (j = 0; j < Ny; j++){
				for (i = 0; i < Nx; i++){
					if(drop[l] || boundary[l]){
						xm = i - 1 + j * Nx + k * Nx * Ny;
						xp = i + 1 + j * Nx + k * Nx * Ny;
						ym = i + (j - 1) * Nx + k * Nx * Ny;
						yp = i + (j + 1) * Nx + k * Nx * Ny;
						zm = i + j * Nx + (k - 1) * Nx * Ny;
						zp = i + j * Nx + (k + 1) * Nx * Ny;
						if(ndrop[xm] || ndrop[xp] || ndrop[ym] || ndrop[yp] || ndrop[zm] || ndrop[zp]){
							if(drop[l]){
								nboundary[l] = true;
								drop[l] = false;
								bulk --;
								nsurf ++;	
								surf ++;
							}	
							else if(Wp > W){
								boundary[l] = false;
								nboundary[l] = true;
								nsurf ++;
							}
						}
					}
					l ++;
				}
			}
		}
	}

	dV = (((4. / 3.) * (double)M_PI * (Rx * Ry * Rz)) - (4.0 / 3.) * (double)M_PI * Rp * Rp * Rp * Np) / bulk; 
	if(DoubleU){
		dVi = (4.0 / 3. * (double)M_PI * (iRx * iRy * iRz)) / (double)bulkin; 
		dVo = (4.0 / 3. * (double)M_PI * ((Rx * Ry * Rz) - (iRx * iRy * iRz))) / (double)(bulkout - surf);
	}

	dAdrop = (4. * (double)M_PI * pow((pow(Rx * Ry, 1.6075) + pow(Rx * Rz, 1.6075) + pow(Ry * Rz, 1.6075)) / 3.0, 1.0/1.6075) - (double)M_PI * Rp * Rp * Np) / (double)(surf -  nsurf);
	if(nsurf > 0){
		dApart = (2 * (double)M_PI * Rp * Rp * Np) / (double)nsurf ;
	}

	else{
		dApart = 0;
	}
//	printf("\ndV is %lf\ndA of droplet is %lf\ndA of nanoparticle is %lf\n", dV, dAdrop, dApart); 

	//calculate dV and dA for NP inside droplet
	//dV = (4.0 / 3 * M_PI * R * R * R - 4.0 / 3 * M_PI * Rp * Rp * Rp * Np) / bulk; 
	//dAdrop = (4 * M_PI * R * R) / (surf -  nsurf);
	//if(Np != 0){
	//	dApart = (4 * M_PI * Rp * Rp * Np) / nsurf ;
	//}
	//else{
	//	dApart = 0;
	//}

	printf("\ndV is %lf\ndA of droplet is %lf\ndA of nanoparticle is %lf\n", dV, dAdrop, dApart); 
	droplet = bulk + surf;
	printf("\nRx is %lf\nRy is %lf\nRz is %lf\nDroplet nodes number is %d\nBulk nodes number is %d\nDroplet surface nodes number is %d\n", Rx, Ry, Rz, droplet, bulk, surf); 
	printf("Internal nodes count = %d\nExternal nodes count = %d\n", bulkin, bulkout);
	printf("External bulk count contains surface nodes!!\n");
	printf("External count after substracting surface nodes is %d\n", bulkout - surf);
	printf("dVi = %lf\ndVo = %lf\n", dVi, dVo);
	//allocate nu 
	//allocate Qo only for finite homeotropic
	nu = (double*)malloc(surf * 3 * sizeof(double));
	for(i = 0; i < surf * 3; i ++){
		nu[i] = 0;
	}
	AnchNInf = false;
	for(int i = 0; i < Np; i++){
		if(pos[i][3] == 1){
			AnchNInf = true;
			break;
		}
	}
	if((degenerate == 0 && infinite == 0) || AnchNInf){
		Qo = (double*)malloc(6 * surf * sizeof(double));
		for(int i = 0; i < surf * 6; i ++){
			Qo[i] = 0;
		}
	}

	//allocate qold and neighbor
	//allocate share to define droplet: -1 not defined; 20 bulk; 0-9 droplet boundary; 10 -19 nanoparticle boundary
	length = lrint(droplet / numprocs) + 1;
	share = (char*)malloc(numprocs * length *  sizeof(char));
	Qold = (double*)malloc(6 * numprocs * length  * sizeof(double));
	neighbor = (int*)malloc(6 * numprocs * length  * sizeof(int));
	for(int i = 0; i < length * numprocs; i ++){
		share[i] = -1;
	}
	for(int i = 0; i < 6 * length * numprocs; i ++){
		Qold[i] = 0;
		neighbor[i] = -1;
	}

	//populate indx array to transformation from 3D to 1D and share array.
	nb = 0;
	nd = 0;
	for(int l = 0; l < tot; l++){
		if(drop[l] || boundary[l] || nboundary[l]){
			indx[l] = nd;
			if(drop[l]) share[nd] = 0;
			else if(boundary[l] || nboundary[l]){
				if(boundary[l])	share[nd] = 2;
				else	share[nd] = 4;
				nb ++;
			}
			nd ++;
		}
	}

	if (nd != droplet){
		printf("Problem in initialization of qtensor. nd is %d not equal to droplet %d.\n", nd, droplet);
		return false;
	}
	if (nb != surf){
		printf("Problem in initialization of qtensor. nb is %d not equal to surf %d.\n", nb, surf);
		return false;
	}

	//Initial configuration of qtensor according to seed 
	if(!conf(pos))	return false;

	//define neighbors and calculate normal vector nu for droplet and particle.
	l = 0;
	nb = 0;
	for(int k = 0; k < Nz; k++){
		for (int j = 0; j < Ny; j++){
			for (int i = 0; i < Nx; i++){
				nd = indx[l];
				if(drop[l]){
					neighbor[nd * 6 + 0] = indx[i - 1 + j * Nx + k * Nx * Ny];
					neighbor[nd * 6 + 1] = indx[i + 1 + j * Nx + k * Nx * Ny];
					neighbor[nd * 6 + 2] = indx[i + (j - 1) * Nx + k * Nx * Ny];
					neighbor[nd * 6 + 3] = indx[i + (j + 1) * Nx + k * Nx * Ny];
					neighbor[nd * 6 + 4] = indx[i + j * Nx + (k - 1) * Nx * Ny];
					neighbor[nd * 6 + 5] = indx[i + j * Nx + (k + 1) * Nx * Ny];
				}
				if(boundary[l] || nboundary[l]){
					if(boundary[l]){
						x = (double)(i - rx) * dx;
						y = (double)(j - ry) * dy;
						z = (double)(k - rz) * dz;
						
						dis = sqrt(x*x+y*y+z*z);
						
						//define nu
						if (dis == 0){
							printf("Error in neighbors on boundary.\n");
							return false;
						}
						
						nu[nb * 3 + 0] = -2. * x / (Rx * Rx);
						nu[nb * 3 + 1] = -2. * y / (Ry * Ry);
						nu[nb * 3 + 2] = -2. * z / (Rz * Rz);
						norm_v(&nu[nb * 3]);
						
						if(DoubleU){

							//infinite, define qtensor and don't evolve any more
							//homeotropic noninfinite, define qo
							if(infinite == 1){
				
								for(n = 0; n < 6; n ++){
									Qold[nd * 6 + n] = dir2ten(&nu[nb * 3], n, S2);
								}	

							}
							else if(degenerate == 0 && infinite == 0){
								for(n = 0; n < 6; n ++){
									Qo[nb * 6 + n] = dir2ten(&nu[nb * 3], n, S2);
								}		
							}

						}
						else{
							//infinite, define qtensor and don't evolve any more
							//homeotropic noninfinite, define qo
							if(infinite == 1){
				
								for(n = 0; n < 6; n ++){
									Qold[nd * 6 + n] = dir2ten(&nu[nb * 3], n, S);
								}	

							}
							else if(degenerate == 0 && infinite == 0){
								for(n = 0; n < 6; n ++){
									Qo[nb * 6 + n] = dir2ten(&nu[nb * 3], n, S);
								}		
							}
						}
					}	
					
					else if(nboundary[l]){
						for(m = 0; m < Np; m ++){
							x = (i- pos[m][0])*dx ;
							y = (j- pos[m][1])*dy ;
							z = (k- pos[m][2])*dz ;
							dis = x*x+y*y+z*z;
							if (dis <= (Rp + 1) * (Rp + 1)){
								if (dis == 0){
									printf("Error in neighbors on particle boundary.\n");
									return false;
								}
								else {
									nu[nb * 3 + 0] = x;
									nu[nb * 3 + 1] = y;
									nu[nb * 3 + 2] = z;
									norm_v(&nu[nb * 3]);
								}

								//PENDIENTE DE REVISAR PARA CONDICIONES DE NANOPARTÍCULAS
								// *************************************************** //
								if(pos[m][3] == 0){
									//Infinite, do not evolve;
									share[nd] = 8;
									for(n = 0; n < 6; n ++){
										Qold[nd * 6 + n] = dir2ten(&nu[nb * 3], n, S); 
									}		
								}
								else if(pos[m][3] == 1){
									//Non-Inf H
									share[nd] = 4;
									for(n = 0; n < 6; n ++){
										Qo[nb * 6 + n] = dir2ten(&nu[nb * 3], n, S);
									}		
								}
								else if(pos[m][3] == 2){
									//Degenerate
									share[nd] = 6;
								}

								// *** PENDIENTE DE REVISAR PARA NANOPARTÍCULAS *** //
								break;
							}
						}

						nsurf --;
					}
					//define boundary
					if(nu[nb * 3 + 0] >= 0){
						neighbor[nd * 6 + 0] = indx[i + 1 + j * Nx + k * Nx * Ny];
						neighbor[nd * 6 + 1] = indx[i + 2 + j * Nx + k * Nx * Ny];
					}
					else if(nu[nb * 3 + 0] < 0){
						neighbor[nd * 6 + 0] = indx[i - 1 + j * Nx + k * Nx * Ny];
						neighbor[nd * 6 + 1] = indx[i - 2 + j * Nx + k * Nx * Ny];
					}
					if(nu[nb * 3 + 1] >= 0){
						neighbor[nd * 6 + 2] = indx[i + (j + 1) * Nx + k * Nx * Ny];
						neighbor[nd * 6 + 3] = indx[i + (j + 2) * Nx + k * Nx * Ny];
					}
					else if(nu[nb * 3 + 1] < 0){
						neighbor[nd * 6 + 2] = indx[i + (j - 1) * Nx + k * Nx * Ny];
						neighbor[nd * 6 + 3] = indx[i + (j - 2) * Nx + k * Nx * Ny];
					}
					if(nu[nb *3 + 2] >= 0){
						neighbor[nd * 6 + 4] = indx[i + j * Nx + (k + 1) * Nx * Ny];
						neighbor[nd * 6 + 5] = indx[i + j * Nx + (k + 2) * Nx * Ny];
					}
					else if(nu[nb * 3 + 2] < 0){
						neighbor[nd * 6 + 4] = indx[i + j * Nx + (k - 1) * Nx * Ny];
						neighbor[nd * 6 + 5] = indx[i + j * Nx + (k - 2) * Nx * Ny];
					}
					nb ++;
				}
				l ++;
			}
		}
	}
	if (nb != surf){
		printf("Problem in initialization of share. nb is %d not equal to surf %d.\n", nb, surf);
		return false;
	}

	for(nd = 0; nd < droplet; nd ++){
		//for all Bulk point, if one of the neighbor is surface point
		count1 = 0;
		if(share[nd] == 0){
			for(n = 0; n < 6; n ++){
				if(share[neighbor[nd * 6 + n]] >= 2){
					count1 ++;
				}
			}
			if(count1 > 1){
				share[nd] += 1;
			} 
		}
		//for all surface point, if one of the neighbor is not defined
		else if(share[nd] < 8 && share[nd] >= 2){	
			for(n = 0; n < 6; n++){
				if(neighbor[nd * 6 + n] == -1){
					count1 ++;	
				}
			}
			if(count1 > 0){
				share[nd] += 1;
			} 
		}
		//for all nodes with problem, share +1
	}

	if(DoubleU){
		int btcount = 0;
		bulktype = (int*)malloc(length * numprocs * sizeof(int));
		for(int i = 0; i < length * numprocs; i++){
			bulktype[i]=0;
			btcount++;
		}

		int bt1 = 0, bt2 = 0;
		nd = 0;
		for(int i = 0; i < tot; i++){
			if(init_bulktype[i] == 1 ){
				bulktype[nd] = 1;
				bt1++;
				nd++;
			}
			else if(init_bulktype[i] == 2){
				bulktype[nd] = 2;
				bt2++;
				nd++;
			}
		}

		if ((bt1 + bt2) != droplet || btcount != length * numprocs) {
			printf("Error in transfer data to droplet bulktype!\n");
			printf("count is %d and droplet size is %d\n", btcount, droplet);
			return false;
		}
		else{
			printf("Data transfer to droplet bulktype successfully!\n");
		}
	}

	free(init_bulktype);
	free(ndrop);
	free(indx);

	if(Np != 0){
		for(m = 0; m < Np; m ++)	free(pos[m]);
		free(pos);
	}

	printf("Initialization of ellipse successful.\n");
	return true;
}
