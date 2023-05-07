#include"finite.h"

int pRx, pRy, pRz;
double pU;
double alpha, beta, gama;
int interface;
int anchoring;

bool read_nano(){

	FILE* param = fopen("nano.in", "r");

	if (param == (FILE*)NULL) {
		printf("No nano.in file found!\n");
		return false;
	}

    fscanf(param, "pRx %d #size of nanoparticle\n", &pRx);
    fscanf(param, "pRy %d\n", &pRy);
    fscanf(param, "pRz %d\n", &pRz);
    fscanf(param, "pU %lf #U for the interface layer\n", &pU);
    fscanf(param, "alpha %lf #Angles for rotate or tilt the nano particle\n", &alpha);
    fscanf(param, "beta %lf\n", &beta);
    fscanf(param, "gamma %lf\n", &gama);
    fscanf(param, "interface %d #thickness of the interface layer; 0: no interface\n", &interface);
    fscanf(param, "anchoring %d #0:random 1:homeotropic 2:planar\n", &anchoring);

    printf("\n~ Nanoparticle data ~\n");
    printf("pRx %d\n", pRx);
    printf("pRy %d\n", pRy);
    printf("pRz %d\n", pRz);
    printf("pU %lf\n", pU);
    printf("alpha: %lf beta: %lf gamma: %lf\n", alpha, beta, gama);
    printf("interface nodes %d\n", interface);
    if(anchoring == 0){
        printf("random anchoring\n");
    }
    else if(anchoring == 1){
        printf("homeotropic anchoring\n");
    }
    else if(anchoring == 2){
        printf("planar anchoring\n");
    }
    else{
        printf("unknonw anchoring\n");
        return false;
    }


    return true;

}

bool initial_nano_channel(){

    bool *ndrop;
    int *indx;
    int nsurf;
    int l;

    int rx, ry, rz;
    double Rx, Ry, Rz;

    //Mitad de la caja
    rx = lrint(Nx / 2) - 1;
    ry = lrint(Ny / 2) - 1;
    rz = lrint(Nz / 2) - 1;

    //Radio del sistema
    Rx = Lx / 2. - 2.;
    Ry = Ly / 2. - 2.;
    Rz = Lz / 2. - 2.;

    double dx = Lx/(Nx-1);
	double dy = Ly/(Ny-1);
	double dz = Lz/(Nz-1);

	idx = 1 / dx;
	idy = 1 / dy;
	idz = 1 / dz;
	iddx = idx * idx;
	iddy = idy * idy;
	iddz = idz * idz;

    surf = 2 * Nx * Ny;
	nsurf = 0;
	tot = Nx * Ny * Nz;
	bulk = Nx * Ny * (Nz - 2);

    //allocate drop and boundary
	ndrop = (bool*)malloc(tot * sizeof(bool));
	nboundary = (bool*)malloc(tot * sizeof(bool));
	drop = (bool*)malloc(tot * sizeof(bool));
	boundary = (bool*)malloc(tot * sizeof(bool));
	indx = (int*)malloc(tot * sizeof(int));
	init_bulktype = (int*)malloc(tot * sizeof(int));
	for(l = 0; l < tot; l ++){
		drop[l] = true;
		boundary[l] = false;
		ndrop[l] = false;
		nboundary[l] = false;
		indx[l] = -1;
	}

    //define the channel surface 
	for (int j = 0; j < Ny; j++){
		for (int i = 0; i < Nx; i++){
			for(int k = 0; k <= Nz - 1; k += Nz - 1){
				l = i + j * Nx + k * Nx * Ny;
				drop[l] = false;
				boundary[l] = true;
			}
		}
	}

    //Reading Nano.in file.
    if(!read_nano()){
        return false;
    }

    alpha = (alpha * M_PI) / 180.;
    beta = (beta * M_PI) / 180.;
    gama = (gama * M_PI) / 180.;

    /*
    Nano particle is in the center of the channel.
    'til now, there's no need of defining a position.
    */

    double x, y, z;
    double x_rot, y_rot, z_rot;
    double distance;
    l = 0;

    for(int k = 0; k < Nz; k++){
        for(int j = 0; j < Ny; j++){
            for(int i = 0; i < Nx; i++){

                x = i - rx;
                y = j - ry;
                z = k - rz;

                x_rot = x * cos(alpha) * cos(beta) + y * (cos(alpha) * sin(beta) * sin(gama) - sin(alpha) * cos(gama))\
					+ z * (cos(alpha) * sin(beta) * cos(gama) + sin(alpha) * sin(gama));
				y_rot = x * sin(alpha) * cos(beta) + y * (sin(alpha) * sin(beta) *sin(gama) + cos(alpha) * cos(gama))\
					+ z * (sin(alpha) * sin(beta) * cos(gama) - cos(alpha) * sin(gama));
				z_rot = x * -sin(beta) + y * cos(beta) * sin(gama) + z * cos(beta) * cos(gama);

				x = x_rot;
				y = y_rot;
				z = z_rot;

                distance = (x * x) / ((pRx + 0.5) * (pRx + 0.5))\
                    + (y * y) / ((pRy + 0.5) * (pRy + 0.5))\
                    + (z * z) / ((pRz + 0.5) * (pRz + 0.5));
                
                if(distance <= 1){
                    ndrop[l] = true;
                    drop[l] = false;
                    bulk--;
                }

                l++;
            
            }
        }
    }

    //Defining Nano boundaries
    /*
    La superficie de la nanopartícula tiene un anclaje infinito.
    Pero la superficie del canal no necesariamente debe tener un
    anclaje infinito también.
    */
    int xm, xp, ym, yp, zm, zp;
    l = 0;
    for(int k = 0; k < Nz; k++){
        for(int j = 0; j < Ny; j++){
            for(int i = 0; i < Nx; i++){

                if(drop[l]){
                    xm = peri(i - 1, 0) + j * Nx + k * Nx * Ny;
                    xp = peri(i + 1, 0) + j * Nx + k * Nx * Ny;
                    ym = i + peri(j - 1, 1) * Nx + k * Nx * Ny;
                    yp = i + peri(j + 1, 1) * Nx + k * Nx * Ny;
                    zm = i + j * Nx + (k - 1) * Nx * Ny;
                    zp = i + j * Nx + (k + 1) * Nx * Ny;
                    
                    if(ndrop[xm] || ndrop[xp] || ndrop[ym] || ndrop[yp] || ndrop[zm] || ndrop[zp]){
						nboundary[l] = true;
						drop[l] = false;
						bulk--;
						nsurf++;	
						surf++;
					}
                }
                l++;          
            }
        }
    }

    /*
    Definiendo bordes de cristal líquido que no evolucionarán.
    Todo dependerá del grosor de la capa (variable interface)
    */

    l = 0;
    if(interface != 0){
        for(int node = 0; node < interface; node++){
            for(int k = 0; k < Nz; k++){
                for(int j = 0; j < Ny; j++){
                    for(int i = 0; i < Nx; i++){
                        if(drop[l]){
                        xm = peri(i - (node + 1), 0) + j * Nx + k * Nx * Ny;
                        xp = peri(i + (node + 1), 0) + j * Nx + k * Nx * Ny;
                        ym = i + peri(j - (node + 1), 1) * Nx + k * Nx * Ny;
                        yp = i + peri(j + (node + 1), 1) * Nx + k * Nx * Ny;
                        zm = i + j * Nx + (k - (node + 1)) * Nx * Ny;
                        zp = i + j * Nx + (k + (node + 1)) * Nx * Ny;
                        
                        if(nboundary[xm] || nboundary[xp] || nboundary[ym]
                            || nboundary[yp] || nboundary[zm] || nboundary[zp]){
                            drop[l] = false;
                            bulk--;
                            surf++;
                            nsurf++;
                        }
                    }
                    l++;
                    }
                }
            }
        }
    }
    
    dV = (Lx * Ly * Lz - 4. / 3. * M_PI * pRx * pRy * pRz) / bulk;
    dAdrop = (2 * Lx * Ly) / (surf - nsurf);
    dApart = 4. * M_PI * pow((pow(pRx * pRy, 1.6075) + pow(pRx * pRz, 1.6075) + pow(pRy * pRz, 1.6075)) / 3.0, 1.0/1.6075) / nsurf;

    int dAinterface;

    printf("\ndV is %lf\ndA of droplet is %lf\ndA of nanoparticle is %lf\n", dV, dAdrop, dApart); 
    droplet = bulk + surf;
	printf("\nDroplet nodes number is %d.\nBulk nodes number is %d.\nDroplet surface nodes number is %d. \nParticle surface nodes number is %d.\n", droplet, bulk, surf, nsurf); 
    
    nu = (double*)malloc(surf * 3 * sizeof(double));
	for(int i = 0; i < surf * 3; i ++){
		nu[i] = 0;
	}

    if(degenerate == 0 && infinite == 0){
		Qo = (double*)malloc(6 * surf * sizeof(double));
		for(int i = 0; i < surf * 6; i ++){
			Qo[i] = 0;
		}
	}

    length = lrint(droplet / numprocs) + 1;
	share = (char*)malloc(numprocs * length * sizeof(char));
	Qold = (double*)malloc(6 * numprocs * length * sizeof(double));
	neighbor = (int*)malloc(6 * numprocs * length * sizeof(int));
	for(int i = 0; i < numprocs * length; i ++){
		share[i] = -1;
	}
	for(int i = 0; i < 6 * numprocs * length; i ++){
		Qold[i] = 0;
		neighbor[i] = -1;
	}

    int num_node = 0;
    int num_boundary = 0;

    for(int i = 0; i < tot; i++){
        if(!ndrop[i]){
            indx[i] = num_node;
            //bulto: share/sign = 0
            //superficie del canal: share/sign = 2
            //superficie de la nanopartícula: share/sign = 4
            if(drop[i]){
                share[num_node] = 0;
            }
            else if(boundary[i]){
                share[num_node] = 2;
                num_boundary++;
            }
            else if(nboundary[i]){
                share[num_node] = 4;
                num_boundary++;
            }
            num_node++;
        }

    }

    if (num_node != droplet){
		printf("Problem in initialization of qtensor. nd is %d not equal to droplet %d.\n", num_node, droplet);
		return false;
	}
	if (num_boundary != surf){
		printf("Problem in initialization of qtensor. nb is %d not equal to surf %d.\n", num_boundary, surf);
		return false;
	}

    double** pos;
    if(!conf(pos)){
        return false;
    }

    int nd = 0;
    int nb = 0;
    time_t t;
    srand((unsigned) time(&t));
    //Defininiendo los vecinos
    for(int k = 0; k < Nz; k++){
        for(int j = 0; j < Ny; j++){
            for(int i = 0; i < Nx; i++){

                nd = indx[i + j * Nx + k * Nx * Ny];

                if(drop[i + j * Nx + k * Nx * Ny]){
                    neighbor[nd * 6 + 0] = indx[peri(i - 1, 0) + j * Nx + k * Nx * Ny];
					neighbor[nd * 6 + 1] = indx[peri(i + 1, 0) + j * Nx + k * Nx * Ny];
					neighbor[nd * 6 + 2] = indx[i + peri(j - 1, 1) * Nx + k * Nx * Ny];
					neighbor[nd * 6 + 3] = indx[i + peri(j + 1, 1) * Nx + k * Nx * Ny];
					neighbor[nd * 6 + 4] = indx[i + j * Nx + (k - 1) * Nx * Ny];
					neighbor[nd * 6 + 5] = indx[i + j * Nx + (k + 1) * Nx * Ny];
                }
                else if(boundary[i + j * Nx + k * Nx * Ny] || nboundary[i + j * Nx + k * Nx * Ny]){
                    if(boundary[i + j * Nx + k * Nx * Ny]){
                        if(k == 0){
                            nu[nb * 3 + 0] = dir1[0];
                            nu[nb * 3 + 1] = dir1[1];
                            nu[nb * 3 + 2] = dir1[2];
                        }
                        else if(k == Nz - 1){
                            nu[nb * 3 + 0] = dir1[0];
                            nu[nb * 3 + 1] = dir1[1];
                            nu[nb * 3 + 2] = -dir1[2];
                        }
                        else{
                            printf("Error in channel surface.\n");
                            return false;
                        }
                        
                        if(infinite == 1){
                            for(int n = 0; n < 6; n ++){
                                Qold[nd * 6 + n] = dir2ten(&nu[nb * 3], n, S);
                            }
                        }
                        else if(degenerate == 0 && infinite == 0){
                            for(int n = 0; n < 6; n ++){
                                Qo[nb * 6 + n] = dir2ten(&nu[nb * 3], n, S);
                            }		
                        }
                    }
                    else if(nboundary[i + j * Nx + k * Nx * Ny]){
                        x = i - rx;
                        y = j - ry;
                        z = k - rz;

                        x_rot = x * cos(alpha) * cos(beta) + y * (cos(alpha) * sin(beta) * sin(gama) - sin(alpha) * cos(gama))\
                            + z * (cos(alpha) * sin(beta) * cos(gama) + sin(alpha) * sin(gama));
                        y_rot = x * sin(alpha) * cos(beta) + y * (sin(alpha) * sin(beta) *sin(gama) + cos(alpha) * cos(gama))\
                            + z * (sin(alpha) * sin(beta) * cos(gama) - cos(alpha) * sin(gama));
                        z_rot = x * -sin(beta) + y * cos(beta) * sin(gama) + z * cos(beta) * cos(gama);

                        x = x_rot;
                        y = y_rot;
                        z = z_rot;

                        distance = (x * x) / ((pRx + 0.5) * (pRx + 0.5))\
                            + (y * y) / ((pRy + 0.5) * (pRy + 0.5))\
                            + (z * z) / ((pRz + 0.5) * (pRz + 0.5));

                        if (distance == 0){
                            printf("Error in neighbors on particle boundary.\n");
                            return false;
                        }
                        else {

                            if(anchoring == 0){
                                nu[nb * 3 + 0] = (rand() % pRx + 1);
						        nu[nb * 3 + 1] = (rand() % pRy + 1);
						        nu[nb * 3 + 2] = (rand() % pRz + 1);
                                norm_v(&nu[nb * 3]);

                            }
                            else if(anchoring == 1){
                                nu[nb * 3 + 0] = 2. * x / (pRx * pRx);
						        nu[nb * 3 + 1] = 2. * y / (pRy * pRy);
						        nu[nb * 3 + 2] = 2. * z / (pRz * pRz);
						        norm_v(&nu[nb * 3]);
                            }
                            else if(anchoring == 2){
                                printf("Not available yet\n");
                                exit(1);
                            }
                            //La superficie no evoluciona
                            share[nd] = 8;
                            Qold[nd * 6 + 0] = dir2ten(&nu[nb * 3], 0, S);
                            Qold[nd * 6 + 1] = dir2ten(&nu[nb * 3], 1, S);
                            Qold[nd * 6 + 2] = dir2ten(&nu[nb * 3], 2, S);
                            Qold[nd * 6 + 3] = dir2ten(&nu[nb * 3], 3, S);
                            Qold[nd * 6 + 4] = dir2ten(&nu[nb * 3], 4, S);
                            Qold[nd * 6 + 5] = dir2ten(&nu[nb * 3], 5, S);
                            
                        }
                        nsurf--;
                    }

                    if(nu[nb * 3 + 0] >= 0){
                        neighbor[nd * 6 + 0] = indx[peri(i + 1, 0) + j * Nx + k * Nx * Ny];
                        neighbor[nd * 6 + 1] = indx[peri(i + 2, 0) + j * Nx + k * Nx * Ny];
                    }
                    else if(nu[nb * 3 + 0] < 0){
                        neighbor[nd * 6 + 0] = indx[peri(i - 1, 0) + j * Nx + k * Nx * Ny];
                        neighbor[nd * 6 + 1] = indx[peri(i - 2, 0) + j * Nx + k * Nx * Ny];
                    }
                    if(nu[nb * 3 + 1] >= 0){
                        neighbor[nd * 6 + 2] = indx[i + peri(j + 1, 1) * Nx + k * Nx * Ny];
                        neighbor[nd * 6 + 3] = indx[i + peri(j + 2, 1) * Nx + k * Nx * Ny];
                    }
                    else if(nu[nb * 3 + 1] < 0){
                        neighbor[nd * 6 + 2] = indx[i + peri(j - 1, 1) * Nx + k * Nx * Ny];
                        neighbor[nd * 6 + 3] = indx[i + peri(j - 2, 1) * Nx + k * Nx * Ny];
                    }
                    if(nu[nb *3 + 2] >= 0){
                        neighbor[nd * 6 + 4] = indx[i + j * Nx + (k + 1) * Nx * Ny];
                        neighbor[nd * 6 + 5] = indx[i + j * Nx + (k + 2) * Nx * Ny];
                    }
                    else if(nu[nb * 3 + 2] < 0){
                        neighbor[nd * 6 + 4] = indx[i + j * Nx + (k - 1) * Nx * Ny];
                        neighbor[nd * 6 + 5] = indx[i + j * Nx + (k - 2) * Nx * Ny];
                    }
                    nb++;
                }
            }
        }
    }

    if (nb != surf){
		printf("Problem in initialization of share. nb is %d not equal to surf %d.\n", nb, surf);
		return false;
	}
    int count1;

    for(nd = 0; nd < droplet; nd ++){
		//for all Bulk point, if one of the neighbor is surface point
		count1 = 0;
		if(share[nd] == 0){
			for(int n = 0; n < 6; n ++){
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
			for(int n = 0; n < 6; n++){
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

	free(ndrop);
	free(indx);

	printf("Initialization of ellipse successful.\n");
	return true;
}