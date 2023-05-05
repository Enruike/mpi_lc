#include"finite.h"

int pRx, pRy, pRz;
double pU;
double alpha, beta, gamma;
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
    fscanf(param, "gamma %lf\n", &gamma);
    fscanf(param, "interface %d #thickness of the interface layer; 0: no interface\n", &interface);
    fscanf(param, "anchoring %d #0:random 1:homeotropic 2:planar\n", &anchoring);

    printf("\n~ Nanoparticle data ~\n");
    printf("pRx %d\n", pRx);
    printf("pRy %d\n", pRy);
    printf("pRz %d\n", pRz);
    printf("U %d\n", pU);
    printf("alpha: %lf beta: %lf gamma: %lf\n", alpha, beta, gamma);
    printf("interface nodes %d\n", interface);
    if(anchoring == 0){
        printf("random anchoring\n");
    }
    else if(anchoring == 1{
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

initial_nano_channel(){

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

    /*
    Nano particle is in the center of the channel.
    'til now, there's no need of defining a position.
    */

    double x, y, z;
    double x_rot, y_rot, z_rot;
    double distance;
    int l = 0;

    for(int k = 0; k < Nz; k++){
        for(int j = 0; j < Ny; j++){
            for(int i = 0; i < Nx; i++){

                x = i - rx + pRx;
                y = j - ry + pRy;
                z = k - rz + pRz;

                x_rot = x * cos(alpha) * cos(beta) + y * (cos(alpha) * sin(beta) * sin(gamma) - sin(alpha) * cos(gamma))\
					+ z * (cos(alpha) * sin(beta) * cos(gamma) + sin(alpha) * sin(gamma));
				y_rot = x * sin(alpha) * cos(beta) + y * (sin(alpha) * sin(beta) *sin(gamma) + cos(alpha) * cos(gamma))\
					+ z * (sin(alpha) * sin(b) * cos(gamma) - cos(alpha) * sin(gamma));
				z_rot = x * -sin(beta) + y * cos(beta) * sin(gamma) + z * cos(beta) * cos(gamma);

				x = x_rot;
				y = y_rot;
				z = z_rot;

                distance = (x * x) / ((pRx + 0.5) * (pRx + 0.5))\
                    + (y * y) / ((pRy + 0.5) * (pRy + 0.5))\
                    + (z * z) / ((pRz + 0.5) * (pRz + 0.5));
                
                if(dis <= 1){
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
    dApart = (4 * M_PI * pRx * pRy * pRz) / nsurf;

    int dAinterface;

    printf("\ndV is %lf\ndA of droplet is %lf\ndA of nanoparticle is %lf\n", dV, dAdrop, dApart); 
    droplet = bulk + surf;
	printf("\nDroplet nodes number is %d.\nBulk nodes number is %d.\nDroplet surface nodes number is %d. \nParticle surface nodes number is %d.\n", droplet, bulk, surf, nsurf); 
    
    nu = (double*)malloc(surf * 3 * sizeof(double));
	for(i = 0; i < surf * 3; i ++){
		nu[i] = 0;
	}
    
    return true;
}