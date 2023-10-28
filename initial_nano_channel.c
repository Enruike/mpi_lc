#include"finite.h"

double pU;
double alpha, beta, gama;
int posX, posY, posZ;

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
    fscanf(param, "posX %d #0 for center; nanoparticle position\n", &posX);
    fscanf(param, "posY %d\n", &posY);
    fscanf(param, "posZ %d\n", &posZ);
    fscanf(param, "pivot %d\n #0:center; 1:edge", &pivotflag);

    if(pivotflag == 0 && (posX != 0 || posY != 0 || posZ != 0)){
        printf("Pivot flag it's set up for 0:center\n");
        printf("Nanoparticle position will be set to center\n");
        posX = 0;
        posY = 0;
        posZ = 0;
    }

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
    if(pivotflag == 0){
        printf("Nanoparticle position will be in the center of the box\n");
    }
    else{
        printf("Nanoparticle position is: posX: %d; posY: %d; posZ: %d\n", posX, posY, posZ);
    }

    return true;

}

bool initial_nano_channel(){

    //Reading Nano.in file.
    if(!read_nano()){
        return false;
        exit(1);
    }

    bool *ndrop;
    int *indx;
    int l;

    //Mitad de la caja
    if(posX == 0){
        rx = lrint(Nx / 2);
    }
    else{
        rx = posX;
    }
    if(posY == 0){
        ry = lrint(Ny / 2);
    }
    else{
        ry = posY;
    }
    if(posZ == 0){
        rz = lrint(Nz / 2);
    }
    else{
        rz = posZ;
    }

    

    //Radio del sistema
    double Rx = Lx / 2. - 2.;
    double Ry = Ly / 2. - 2.;
    double Rz = Lz / 2. - 2.;

    dx = Lx/(Nx-1);
	dy = Ly/(Ny-1);
	dz = Lz/(Nz-1);

	idx = 1 / dx;
	idy = 1 / dy;
	idz = 1 / dz;
	iddx = idx * idx;
	iddy = idy * idy;
	iddz = idz * idz;

    surf = 0;
	nsurf = 0;
	tot = Nx * Ny * Nz;
	bulk = Nx * Ny * Nz;

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
        init_bulktype[l] = 1;
		indx[l] = 1;
	}

 
    //define the channel surface 
	for(int k = 0; k < Nz; k += (Nz - 1)){
        for(int j = 0; j < Ny; j++){
            for(int i = 0; i < Nx; i++){
                l = i + j * Nx + k * Nx * Ny;
                boundary[l] = true;
                drop[l] = false;
                init_bulktype[l] = 4;
                surf++;
                bulk--;
            }
        }
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
    int nanoparticle_nodes = 0;
    l = 0;

    double pivotX;
    double pivotY;
    double pivotZ;

    //Pivot point.
    //pivotX = (double)pRx;
    //pivotY = (double)posY;
    //pivotZ = (double)posZ;
    pivotX = 0.;
    pivotY = 0.;
    pivotZ = sin(beta) * pRx;

    for(int k = 0; k < Nz; k++){
        for(int j = 0; j < Ny; j++){
            for(int i = 0; i < Nx; i++){

                x = (double)(i - rx) * dx;
				y = (double)(j - ry) * dy;
				z = (double)(k - rz) * dz;

                //pivotflag 0 for center
                if(pivotflag == 0){
                    x_rot = x * cos(alpha) * cos(beta) + y * (cos(alpha) * sin(beta) * sin(gama) - sin(alpha) * cos(gama))\
					    + z * (cos(alpha) * sin(beta) * cos(gama) + sin(alpha) * sin(gama));
				    y_rot = x * sin(alpha) * cos(beta) + y * (sin(alpha) * sin(beta) *sin(gama) + cos(alpha) * cos(gama))\
					    + z * (sin(alpha) * sin(beta) * cos(gama) - cos(alpha) * sin(gama));
				    z_rot = x * -sin(beta) + y * cos(beta) * sin(gama) + z * cos(beta) * cos(gama);
                }
                else{
                    x_rot = (x - pivotX) * cos(alpha) * cos(beta) + (y - pivotY) * (cos(alpha) * sin(beta) * sin(gama) - sin(alpha) * cos(gama))\
					    + (z - pivotZ) * (cos(alpha) * sin(beta) * cos(gama) + sin(alpha) * sin(gama));
				    y_rot = (x - pivotX) * sin(alpha) * cos(beta) + (y - pivotY)* (sin(alpha) * sin(beta) *sin(gama) + cos(alpha) * cos(gama))\
					    + (z - pivotZ) * (sin(alpha) * sin(beta) * cos(gama) - cos(alpha) * sin(gama));
				    z_rot = (x - pivotX) * -sin(beta) + (y - pivotY) * cos(beta) * sin(gama) + (z - pivotZ) * cos(beta) * cos(gama);
                }

				x = x_rot;
				y = y_rot;
				z = z_rot;

                distance = (x * x) / ((pRx + 0.5) * (pRx + 0.5))\
                    + (y * y) / ((pRy + 0.5) * (pRy + 0.5))\
                    + (z * z) / ((pRz + 0.5) * (pRz + 0.5));
                
                if(distance <= 1){
                    ndrop[l] = true;
                    nanoparticle_nodes++;
                    drop[l] = false;
                    init_bulktype[l] = 5;
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
                        init_bulktype[l] = 6;
						bulk--;
						nsurf++;
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
    int interbulk = 0;

    if(interface != 0){
       
        int rebulker = bulk;

        
        for(int k = 0; k < Nz; k++){
            for(int j = 0; j < Ny; j++){
                for(int i = 0; i < Nx; i++){
                    if(nboundary[l]){

                        for(int node = 0; node < interface; node++){
                            xm = peri(i - (node + 1), 0) + j * Nx + k * Nx * Ny;
                            xp = peri(i + (node + 1), 0) + j * Nx + k * Nx * Ny;
                            ym = i + peri(j - (node + 1), 1) * Nx + k * Nx * Ny;
                            yp = i + peri(j + (node + 1), 1) * Nx + k * Nx * Ny;
                            zm = i + j * Nx + peri(k - (node + 1), 2) * Nx * Ny;
                            zp = i + j * Nx + peri(k + (node + 1), 2) * Nx * Ny;
                            
                            
                            if(drop[xm] && init_bulktype[xm] != 3){
                                init_bulktype[xm] = 3;
                                interbulk++;
                                rebulker--;
                            }
                            if(drop[xp] && init_bulktype[xp] != 3){
                                init_bulktype[xp] = 3;
                                interbulk++;
                                rebulker--;
                            }
                            if(drop[ym] && init_bulktype[ym] != 3){
                                init_bulktype[ym] = 3;
                                interbulk++;
                                rebulker--;
                            }
                            if(drop[yp] && init_bulktype[yp] != 3){
                                init_bulktype[yp] = 3;
                                interbulk++;
                                rebulker--;
                            }
                            if(drop[zm] && init_bulktype[zm] != 3){
                                init_bulktype[zm] = 3;
                                interbulk++;
                                rebulker--;
                            }
                            if(drop[zp] && init_bulktype[zp] != 3){
                                init_bulktype[zp] = 3;
                                interbulk++;
                                rebulker--;
                            }
                        } 
                    }
                l++;
                }
            }
        }
        

        if(bulk != (rebulker + interbulk)){
            printf("Problems with interface nodes\n");
            return false;
        }
    }

    
    
    dV = (Lx * Ly * Lz - 4. / 3. * M_PI * pRx * pRy * pRz) / bulk;
    dVi = (Lx * Ly * Lz - 4. / 3. * M_PI * pRx * pRy * pRz) / (bulk - interbulk);
    if(interface != 0) dVo = (4. / 3. * M_PI * ((double)(pRx + interface) * (double)(pRy + interface) * (double)(pRz + interface) - (double)(pRx) * (double)(pRy) * (double)(pRz))) / (double)interbulk;
    else dVo = 0.;
    dAdrop = (2 * Lx * Ly) / (surf);
    dApart = 4. * M_PI * pow((pow(pRx * pRy, 1.6075) + pow(pRx * pRz, 1.6075) + pow(pRy * pRz, 1.6075)) / 3.0, 1.0/1.6075) / nsurf;

    int dAinterface;

    droplet = bulk + surf + nsurf;

    printf("\ndV is %lf\ndA of droplet is %lf\ndA of nanoparticle is %lf\n", dV, dAdrop, dApart); 
    printf("dVi = %lf\ndVo = %lf\n", dVi, dVo);
	printf("\nDroplet nodes number is %d.\nBulk nodes number is %d.\nDroplet surface nodes number is %d.\nNanoparticle nodes is %d.\nParticle surface nodes number is %d.\n", droplet, bulk, surf, nanoparticle_nodes, nsurf);
    printf("Nanoparticle interface is %d\n", interbulk);

    nu = (double*)malloc((surf + nsurf) * 3 * sizeof(double));
	for(int i = 0; i < (surf + nsurf) * 3; i ++){
		nu[i] = 0;
	}

    if(degenerate == 0 && infinite == 0){
		Qo = (double*)malloc(6 * (surf + nsurf) * sizeof(double));
		for(int i = 0; i < (surf + nsurf) * 6; i ++){
			Qo[i] = 0;
		}
	}

    length = lrint(droplet / numprocs) + 1;
	share = (signed char*)malloc(numprocs * length * sizeof(signed char));
	Qold = (double*)malloc(6 * numprocs * length * sizeof(double));
	neighbor = (int*)malloc(6 * numprocs * length * sizeof(int));
	for(int i = 0; i < numprocs * length; i ++){
		share[i] = -1;
	}
	for(int i = 0; i < 6 * numprocs * length; i ++){
		Qold[i] = 0;
		neighbor[i] = -1;
	}

    int nd = 0;
    int nb = 0;
    int nbulk = 0;

    for(int i = 0; i < tot; i++){
        if(!ndrop[i]){
            indx[i] = nd;
            //bulto: share/sign = 0
            //superficie del canal: share/sign = 2
            //superficie de la nanopartícula: share/sign = 4
            if(drop[i]){
                share[nd] = 0;
                nbulk++;
            }
            else if(boundary[i]){
                share[nd] = 2;
                nb++;
            }
            else if(nboundary[i]){
                share[nd] = 4;
                nb++;
            }
            nd++;
        }

    }
    int count1;
    int countshare0 = 0;
    int countshare2 = 0;
    int countshare4 = 0;
    int countshare8 = 0;
    int shareminusone = 0;
    int undefined = 0;
    count1 = 0;
    for(int i = 0; i < droplet; i++){
        if(share[i] == 0){
            countshare0++;
            count1++;
        }
        else if(share[i] == 2){
            countshare2++;
            count1++;
        }
        else if(share[i] == 4){
            countshare4++;
            count1++;
        }
        else if(share[i] == 8){
            countshare8++;
            count1++;
        }
        else if(share[i] == -1){
            shareminusone++;
        }
        else{
            undefined++;
            
        }
    }
    printf("Pre share count 0 : %d, 2 : %d, 4 : %d, 8 : %d, total : %d\n", countshare0, countshare2, countshare4, countshare8, count1);
    printf("-1 : %d, undefined : %d\n", shareminusone, undefined);

    //printf("nbulk count : %d\n", nbulk);
    if (nd != droplet){
		printf("Problem in initialization of qtensor. nd is %d not equal to droplet %d.\n", nd, droplet);
		return false;
	}
	if (nb != surf + nsurf){
		printf("Problem in initialization of qtensor. nb is %d not equal to surf %d.\n", nb, surf);
        printf("Channel surface nodes: %d; Nanoparticle surface nodes %d; Nanoboundary counter: %d\n", surf, nsurf, nb);
		return false;
	}

    double** pos;
    if(!conf(pos)){
        return false;
    }

    nd = 0;
    nb = 0;
    // time_t t;
    // srand((unsigned) time(&t));
    srand(rand_seed);
    //Defininiendo los vecinos
    for(int k = 0; k < Nz; k++){
        for(int j = 0; j < Ny; j++){
            for(int i = 0; i < Nx; i++){

                nd = indx[i + j * Nx + k * Nx * Ny];

                if(nd == -1){
                    continue;
                }

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
                
                        norm_v(&nu[nb * 3]);

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
                        x = (double)(i - rx) * dx;
                        y = (double)(j - ry) * dy;
                        z = (double)(k - rz) * dz;

                        if(pivotflag == 0){
                            x_rot = x * cos(alpha) * cos(beta) + y * (cos(alpha) * sin(beta) * sin(gama) - sin(alpha) * cos(gama))\
                                + z * (cos(alpha) * sin(beta) * cos(gama) + sin(alpha) * sin(gama));
                            y_rot = x * sin(alpha) * cos(beta) + y * (sin(alpha) * sin(beta) *sin(gama) + cos(alpha) * cos(gama))\
                                + z * (sin(alpha) * sin(beta) * cos(gama) - cos(alpha) * sin(gama));
                            z_rot = x * -sin(beta) + y * cos(beta) * sin(gama) + z * cos(beta) * cos(gama);
                        }
                        else{
                            x_rot = (x - pivotX) * cos(alpha) * cos(beta) + (y - pivotY) * (cos(alpha) * sin(beta) * sin(gama) - sin(alpha) * cos(gama))\
                                + (z - pivotZ) * (cos(alpha) * sin(beta) * cos(gama) + sin(alpha) * sin(gama));
                            y_rot = (x - pivotX) * sin(alpha) * cos(beta) + (y - pivotY)* (sin(alpha) * sin(beta) *sin(gama) + cos(alpha) * cos(gama))\
                                + (z - pivotZ) * (sin(alpha) * sin(beta) * cos(gama) - cos(alpha) * sin(gama));
                            z_rot = (x - pivotX) * -sin(beta) + (y - pivotY) * cos(beta) * sin(gama) + (z - pivotZ) * cos(beta) * cos(gama);
                        }

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
                                nu[nb * 3 + 0] = (rand() % pRx - pRx);
						        nu[nb * 3 + 1] = (rand() % pRy - pRy);
						        nu[nb * 3 + 2] = (rand() % pRz - pRz);
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

    if (nb != surf + nsurf){
		printf("Problem in initialization of share. nb is %d not equal to surf %d.\n", nb, surf);
		return false;
	}
    
    countshare0 = 0;
    countshare2 = 0;
    countshare4 = 0;
    countshare8 = 0;
    shareminusone = 0;
    undefined = 0;
    count1 = 0;
    for(int i = 0; i < droplet; i++){
        if(share[i] == 0){
            countshare0++;
            count1++;
        }
        else if(share[i] == 2){
            countshare2++;
            count1++;
        }
        else if(share[i] == 4){
            countshare4++;
            count1++;
        }
        else if(share[i] == 8){
            countshare8++;
            count1++;
        }
        else if(share[i] == -1){
            shareminusone++;
        }
        else{
            undefined++;
        }
    }

    printf("After count 0 : %d, 2 : %d, 4 : %d, 8 : %d, total : %d\n", countshare0, countshare2, countshare4, countshare8, count1);
    printf("-1 : %d, undefined : %d\n", shareminusone, undefined);

    if(DoubleU){

		bulktype = (int*)malloc(length * numprocs * sizeof(int));
		for(int i = 0; i < length * numprocs; i++){
			bulktype[i]=0;
		}

		int bt1 = 0, bt2 = 0, bt0 = 0;
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
            else if(init_bulktype[i] == 3){
                bulktype[nd] = 3;
                share[nd] = 10; // 10 for interface type. Share window
                bt2++;
                nd++;
               
            }
            else if(init_bulktype[i] == 4 || init_bulktype[i] == 6){
                bulktype[nd] = init_bulktype[i];
                bt2++;
                nd++;
            }
            else if(init_bulktype[i] == -1){
                bt0++;
                nd++;
            } 
		}

		if ((bt1 + bt2) != droplet) {
			printf("Error in transfer data to droplet bulktype!\n");
            printf("bt1 is %d bt2 is %d, bt0 is %d\n", bt1, bt2, bt0);
			printf("droplet size is %d\n", droplet);
			return false;
		}
		else{
			printf("Data transfer to droplet bulktype successfully!\n");
		}
	}

    if (interface != 0)
    {
        int sharecounter = 0;
        for(int i = 0; i < numprocs * length; i++){
            if (share[i] == 10)
            {
                sharecounter++;
            }
            
        }

        if (sharecounter != interbulk)
        {
            printf("Problems in share and interface nodes count!\n");
            printf("Share count is %d and interbulk is %d\n", sharecounter, interbulk);
            exit(1);
        }
        else{
            printf("Share count for interface nodes is ok\n");
        }
        
    }
    
    
    

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

    count1 = 0;

    for(int i = 0; i < droplet; i++){
        if(share[i] >= 0 && share[i] < 14){
            count1++;
        }
    }
    printf("Final count before scattering is %d\n", count1);

	free(ndrop);
	free(indx);
    free(init_bulktype);

	printf("Initialization of nanocanal successful.\n");
	return true;
}