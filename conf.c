#include "finite.h" 

bool conf(double **pos){
	int l;
	int nb, nd;
	int i, j, k, n;
	double dis, x = 0.0, y = 0.0, z = 0.0, xi, yi, zi;
	double disxy;
	double sinthe, costhe, sinphi, cosphi, omega;
	double Qini[6] = {0};
	double dir[3] = {0};
	int flag = 0;
	bool ideal = false;

	//uniform initial configuration
	if(seed == 0){

		printf("Uniform configuration.\n");
		nd = 0;
		double dirvec1[3] = {0};
		double dirvec2[3] = {0};
		double norm = 0;
		srand(rand_seed);
		double dir_temp[3] = { 0. };
		l = 0;
		
		int tempcounter = 0;

		if(geo == -2 || geo == -3){

			for(int k = 0; k < Nz; k++){
				for(int j = 0; j < Ny; j++){
					for(int i = 0; i < Nx; i++){

						if(drop[i + j * Nx + k * Nx * Ny]){						

							for(int n = 0; n < 6; n++){
									Qold[nd * 6 + n] = dir2ten(init_dir, n, S);
							}						
							nd++;
							
						}

						else if(boundary[i + j * Nx + k * Nx * Ny]){
							
							if(j == 2){

								for(int n = 0; n < 6; n++){
									Qold[nd * 6 + n] = dir2ten(dir2, n, S);
								}

							}

							else{

								
								if(uppersurf == 1){

									if(geo == -2){
										x = (i - rx) * dx;
										y = (j - 2) * dy;
										

										norm = sqrt(x * x + y * y);

										x = x/norm;
										y = y/norm;

										dirvec1[0] = -x;
										dirvec1[1] = -y;
										dirvec1[2] = 0;
										
										for(int n = 0; n < 6; n++){
											Qold[nd * 6 + n] = dir2ten(dirvec1, n, S);
										}										
									}

									else if(geo == -3){
										x = (i - rx) * dx;
										y = (j - 2) * dy;
										z = (k - rz) * dz;

										norm = sqrt(x * x + y * y + z * z);

										x = x/norm;
										y = y/norm;
										z = z/norm;

										dirvec1[0] = -x;
										dirvec1[1] = -y;
										dirvec1[2] = -z;
										
										for(int n = 0; n < 6; n++){
											Qold[nd * 6 + n] = dir2ten(dirvec1, n, S);
										}										
									}
									

								}

								else{
									
									for(int n = 0; n < 6; n++){
										Qold[nd * 6 + n] = dir2ten(dir1, n, S);
									}					

								}							

							}
							nd++;
						}
					}
				}
			}		
		}
		else if(geo == 10){
			
			if(!norm_v(init_dir)){
				printf("Problems in initial direction before loop!\n");
				exit(1);
			}
			for(int i = 0; i < 6; i++){
				Qini[i] = dir2ten(init_dir, i, S);
			}
			dirvec1[0] = dirvec2[0] = init_dir[0];
			dirvec1[1] = dirvec2[1] = init_dir[1];
			dirvec1[2] = init_dir[2];
			dirvec2[2] = - init_dir[2];
			norm_v(dirvec1);
			norm_v(dirvec2);
			bool flag = true;

			for(int k = 0; k < Nz; k++){
				for(int j = 0; j < Ny; j++){
					for(int i = 0; i < Nx; i++){
						if(drop[l] || boundary[l] || nboundary[l]){

							
							if(init_bulktype[l] == 3){
								
								do{
									//Modified condition for negative direcctions.
									dir_temp[0] = rand() % (2 * pRx + 2 * interface) - pRx - interface;
									dir_temp[1] = rand() % (2 * pRy + 2 * interface) - pRy - interface;
									dir_temp[2] = rand() % (2 * pRz + 2 * interface) - pRz - interface;
									
									/* original condition for positive numbers 
									
									dir_temp[0] = (rand() % (pRx + interface) + 1);
									dir_temp[1] = (rand() % (pRy + interface) + 1);
									dir_temp[2] = (rand() % (pRz + interface) + 1);
									
									*/
									flag = norm_v(dir_temp);

								}
								while(!flag);

								flag = false;

								/* if(tempcounter < 50) {
									printf("Vector is x:%lf y:%lf z:%lf\n", dir_temp[0], dir_temp[1], dir_temp[2]);
									
								} */
								//printf("before norm dir1 %lf dir2 %lf dir3 %lf\n", dir_temp[0], dir_temp[1], dir_temp[2]);
								if(sqrt(pow(dir_temp[0], 2) + pow(dir_temp[1], 2) + pow(dir_temp[2], 2)) == 0.){
									printf("Problems with random directions!\n");
									printf("Problematic vector is x:%lf y:%lf z:%lf\n", dir_temp[0], dir_temp[1], dir_temp[2]);
									exit(1);
								}
								//printf("after norm dir1 %lf dir2 %lf dir3 %lf\n", dir_temp[0], dir_temp[1], dir_temp[2]);

								
							/* 	if(tempcounter < 50) {
									printf("Norm vector is x:%lf y:%lf z:%lf\n", dir_temp[0], dir_temp[1], dir_temp[2]);
									tempcounter++;
								} */

								Qold[nd * 6 + 0] = dir2ten(dir_temp, 0, S2);
								Qold[nd * 6 + 1] = dir2ten(dir_temp, 1, S2);
								Qold[nd * 6 + 2] = dir2ten(dir_temp, 2, S2);
								Qold[nd * 6 + 3] = dir2ten(dir_temp, 3, S2);
								Qold[nd * 6 + 4] = dir2ten(dir_temp, 4, S2);
								Qold[nd * 6 + 5] = dir2ten(dir_temp, 5, S2);
							}

							else{

								
								if(k == 0){
									for(int m = 0; m < 6; m++){
										Qold[nd * 6 + m] = dir2ten(dirvec1, m, S);
									}
								}
								else if(k == Nz - 1){
									for(int m = 0; m < 6; m++){
										Qold[nd * 6 + m] = dir2ten(dirvec2, m, S);
									}
								}
								else{
									for(int m = 0; m < 6; m++){
										Qold[nd * 6 + m] = Qini[m];
									}									
								}
	
							}
							nd++;
						}
						l++;
					}
				}
			}
		}

		else{

			if(!norm_v(init_dir)){			
				printf("Problems in initial direction! \n");
				return false;
			}			

			for(n = 0; n < 6; n ++){
				Qini[n] = dir2ten(init_dir, n, S);
			}

			for(nd = 0; nd < droplet; nd ++){
				if(drop[nd]){
					for(n = 0; n < 6; n ++){
						Qold[nd * 6 + n] = Qini[n];
					}
				}
			}
		}    	

	}

	else if(seed == 1){
		srand(rand_seed);
		nd = 0;
		double norm = 0.0;
		double dirvec1[3] = {0};

		if((geo == -2 || geo == -3) && ideal == 1){

			for(int k = 0; k < Nz; k++){
				for(int j = 0; j < Ny; j++){
					for(int i = 0; i < Nx; i++){

						if(drop[i + j * Nx + k * Nx * Ny]){	

							for(int n = 0; n < 3; n ++){
								dir[n] = (double)rand() / (double)RAND_MAX * 2 - 1;
							}

							for(int n = 0; n < 6; n++){
								Qold[nd * 6 + n] = dir2ten(dir, n, S);
							}
							nd++;
						}

						else if(boundary[i + j * Nx + k * Nx * Ny]){
							
							if(j == 2){
								
								for(int n = 0; n < 6; n++){
									Qold[nd * 6 + n] = dir2ten(dir2, n, S);
								}				

							}

							else{		

								if(uppersurf == 1){

									if(geo == -2){
										x = (i - rx) * dx;
										y = (j - 2) * dy;

										norm = sqrt(x * x + y * y);

										x = x/norm;
										y = y/norm;

										dirvec1[0] = -x;
										dirvec1[1] = -y;
										dirvec1[2] = 0;

										if(!norm_v(dirvec1)){
											printf("Problems with upper vector initialization\n");
											return false;
										} 
										
										for(int n = 0; n < 6; n++){
											Qold[nd * 6 + n] = dir2ten(dirvec1, n, S);
										}										
									}

									else if(geo == -3){

										x = (i - rx) * dx;
										y = (j - 2) * dy;
										z = (k - rz) * dz;

										norm = sqrt(x * x + y * y + z * z);

										x = x/norm;
										y = y/norm;
										z = z/norm;

										dirvec1[0] = -x;
										dirvec1[1] = -y;
										dirvec1[2] = -z;
										
										if(!norm_v(dirvec1)){
											printf("Problems with upper vector initialization\n");
											return false;
										} 

										for(int n = 0; n < 6; n++){
											Qold[nd * 6 + n] = dir2ten(dirvec1, n, S);
										}										
									}
									

								}

								else{	

									if(!norm_v(dir1)){
										printf("Problems with upper vector initialization\n");
										return false;
									} 	

									for(int n = 0; n < 6; n++){
										Qold[nd * 6 + n] = dir2ten(dir1, n, S);
									}					

								}							

							}
							nd++;
						}
					}
				}
			}		
		}

		else{
			for(int nd = 0; nd < droplet; nd ++){
				for(n = 0; n < 3; n ++){
					dir[n] = (double)rand() / (double)RAND_MAX * 2 - 1;
				}

				if(!norm_v(dir)){
					printf("Problems with random initialization\n");
					return false;
				}        	

				for(n = 0; n < 6; n ++){
					Qold[nd * 6 + n] = dir2ten(dir, n, 0.5);
				}            
			}
		}
	}

	//random near particle
	else if(seed == 10){
		l = 0;
		nd = 0;
		srand(rand_seed);

		if(!norm_v(init_dir))   return false;
			for(k = 0; k < Nz; k++){
				for (j = 0; j < Ny; j++){
					for (i = 0; i < Nx; i++){
						if(drop[l] || boundary[l] || nboundary[l]){
							dir[0] = init_dir[0];			
							dir[1] = init_dir[1];			
							dir[2] = init_dir[2];			
							for(n = 0; n < Np; n++){
								x = i - pos[n][0];
								y = j - pos[n][1];
								z = k - pos[n][2];
								if(x > 0.5 * Nx){
									x -= Nx;
								}
								else if(x < -0.5 * Nx){
									x += Nx;
								}
								if(y > 0.5 * Ny){
									y -= Ny;
								}
								else if(y < -0.5 * Ny){
									y += Ny;
								}
								dis = sqrt(x*x+y*y+z*z);
								if(dis < 2 * Rp){
									for(n = 0; n < 3; n ++){
										dir[n] = (double)rand() / (double)RAND_MAX * 2 - 1;
									}
								}

							}
							if(!norm_v(dir))   return false;
							for(n = 0; n < 6; n ++){
								Qold[nd * 6 + n] = dir2ten(dir, n, S);
							}
							nd ++;
						}
						l ++;
					}
				}
			}
	}

	//DSS or RSS initial configuration
	else if(seed == 2 || seed == 3){
		l = 0;
		nd = 0;
		for(k = 0; k < Nz; k++){
			for (j = 0; j < Ny; j++){
				for (i = 0; i < Nx; i++){
					if(drop[l] || boundary[l]||nboundary[l]){
						
						x = (i-rx)*dx;
						y = (j-ry)*dy;
						z = (k-rz)*dz;
						dis = sqrt(x*x+y*y+z*z);

						if(seed == 2){
							omega = dis * qch;
						}   
						
						else{
							omega = atan2(y, x) + dis * qch;
							disxy = sqrt(x * x + y * y);
						}
							
						if(disxy == 0){
							dir[2] = 1;
							dir[0] = dir[1] = 0;
						}

						else{
							costhe = z / dis;
							sinthe = disxy / dis;
							cosphi = x / disxy;
							sinphi = y / disxy;
							dir[0] = cos(omega) * costhe * cosphi - sin(omega) * sinphi;
							dir[1] = cos(omega) * costhe * sinphi + sin(omega) * cosphi;
							dir[2] = - cos(omega) * sinthe;
						}
		
						if(!norm_v(dir)){
							return false;
						}
							
						for (n = 0; n < 6; n++) {
							Qold[nd * 6 + n] = dir2ten(dir, n, S);
						}
						
						nd ++;
					}
					l ++;
				}
			}
		}
    }

    //initial configuration read from Qtensor_ini.out file
    else if(seed == -1){
        if(!norm_v(init_dir)){
			return false;
		}
		
		for(n = 0; n < 6; n ++){
        	Qini[n] = dir2ten(init_dir, n, 0.5);
        }
                
		double a[6] = {0};
                
		FILE* qtensor;
        qtensor = fopen("Qtensor.bin","rb");
		FILE* grid;
		grid = fopen("grid.bin", "rb");
		int signal;
			
		if(qtensor == (FILE*)NULL){
			printf("File Qtensor.bin not found.\n");
			return false;
		}
			
		if(grid == (FILE*)NULL){
			printf("File grid.bin not found.\n");
			return false;
		}
	
		nd = 0;
		for(l = 0; l < tot; l++){
			fread(&signal, sizeof(int), 1, grid);
			
			if(signal == 0 || signal == 1){
				fread(a, sizeof(double), 6, qtensor);
					a[5] = - a[0] - a[3];
			}				
			else{
				for (n = 0; n < 6; n++) {
					a[n] = Qini[n]; 
				}
			}
			
			if(drop[l] || boundary[l] || nboundary[l]){
				for (n = 0; n < 6; n++) {
					Qold[nd * 6 + n] = a[n];
				}
				nd ++;
			}
			else{
				for(n = 0; n < 6; n ++){
					Qold[nd * 6 + n] = Qini[n];
				}
			}
		}
		
		fclose(qtensor);
		fclose(grid);
}

	else if(seed == -1442 || seed == -1443 || seed == -1444 || seed == -1445 ||
			seed == -1446){
		
		if(!norm_v(init_dir)){
			return false;
		}
		
		for(n = 0; n < 6; n ++){
        	Qini[n] = dir2ten(init_dir, n, 0.5);
        }
                
		double a[6] = {0};
                
		FILE* qtensor;
        qtensor = fopen("Qtensor_shell.bin","rb");
		FILE* grid;
		grid = fopen("grid.bin", "rb");
		int signal;
			
		if(qtensor == (FILE*)NULL){
			printf("File Qtensor_shell.bin not found.\n");
			return false;
			exit(1);
		}
			
		if(grid == (FILE*)NULL){
			printf("File grid.bin not found.\n");
			return false;
		}
	
		nd = 0;

		for(l = 0; l < tot; l++){


			fread(&signal, sizeof(int), 1, grid);
			
			if(signal == 0 || signal == 1){
				fread(a, sizeof(double), 6, qtensor);
					a[5] = - a[0] - a[3];
			}				
			else{
				for (n = 0; n < 6; n++) {
					a[n] = Qini[n]; 
				}
			}

			if(init_bulktype[l] == 1){
				nd++;
			}
			
			else if(init_bulktype[l] == 2 || init_bulktype[l] == 3){
				for (n = 0; n < 6; n++) {
					Qold[nd * 6 + n] = a[n];
				}
				nd ++;
			}
		}
		
		fclose(qtensor);
		fclose(grid);

		l = 0;
		nd = 0;

		if(seed == -1442 || seed == -1443){

			for(int k = 0; k < Nz; k++){
				for(int j = 0; j < Ny; j++){
					for(int i = 0; i < Nx; i++){
						
						if(init_bulktype[l] == 1){
							
							x = (i-rx)*dx;
							y = (j-ry)*dy;
							z = (k-rz)*dz;
							dis = sqrt(x*x+y*y+z*z);

							if(seed == -1442){
								omega = dis * qch;
							}   
							
							else{
								omega = atan2(y, x) + dis * qch;
								disxy = sqrt(x * x + y * y);
							}
								
							if(disxy == 0){
								dir[2] = 1;
								dir[0] = dir[1] = 0;
							}

							else{
								costhe = z / dis;
								sinthe = disxy / dis;
								cosphi = x / disxy;
								sinphi = y / disxy;
								dir[0] = cos(omega) * costhe * cosphi - sin(omega) * sinphi;
								dir[1] = cos(omega) * costhe * sinphi + sin(omega) * cosphi;
								dir[2] = - cos(omega) * sinthe;
							}
			
							if(!norm_v(dir)){
								return false;
							}
								
							for (n = 0; n < 6; n++) {
								Qold[nd * 6 + n] = dir2ten(dir, n, S);
							}
							
							nd++;
						}
						else if(init_bulktype[l] == 2 || init_bulktype[l] == 3){
							nd++;
						}
						l++;
					}
				}
			}
		}
		else if(seed == -1444 || seed == -1445 || seed == -1446){

			double A = 0.2;
			double cst;
			double theta = 45 / 180.0 * M_PI;
			double isq2 = 1.0 / sqrt(2.);
			double sq2 = sqrt(2.);
						
			cst = 2. * qch * redshift;

			if(seed == -1444){

				for(int k = 0; k < Nz; k++){
					for (int j = 0; j < Ny; j++){
						for (int i = 0; i < Nx; i++){
							
							if(init_bulktype[l] == 1){
								
								x = (double)(i - rx) * cst * isq2;
								y = (double)(j - ry) * cst * isq2;
								z = (double)(k - rz) * cst * isq2;
								
											
								Qold[nd * 6 + 0] = A * (- sin(y) * cos(x) - sin(x) * cos(z) + 2 * sin(z) * cos(y));
								Qold[nd * 6 + 3] = A * (- sin(z) * cos(y) - sin(y) * cos(x) + 2 * sin(x) * cos(z));
								Qold[nd * 6 + 5] = A * (- sin(x) * cos(z) - sin(z) * cos(y) + 2 * sin(y) * cos(x));
								Qold[nd * 6 + 1] = A * (- sq2 * sin(x) * sin(z) - sq2 * cos(y) * cos(z) + sin(x) * cos(y));
								Qold[nd * 6 + 2] = A * (- sq2 * sin(z) * sin(y) - sq2 * cos(x) * cos(y) + sin(z) * cos(x));
								Qold[nd * 6 + 4] = A * (- sq2 * sin(y) * sin(x) - sq2 * cos(z) * cos(x) + sin(y) * cos(z));

								nd ++;
							}
							else if(init_bulktype[l] == 2 || init_bulktype[l] == 3){
								nd++;
							}
							l++;
						}
					}
				}
			}

			else if (seed == -1445){

				for(int k = 0; k < Nz; k++){
					for(int j = 0; j < Ny; j++){
						for(int i = 0; i < Nx; i++){
							
							if(init_bulktype[l] == 1){
								
								x = (i - rx) * dx;
								y = (j - ry) * dy;
								z = (k - rz) * dz;
								
								Qold[nd * 6 + 0] = A * (cos(cst * z) - cos(cst * y));
								Qold[nd * 6 + 1] = A * sin(cst * z);
								Qold[nd * 6 + 2] = A * sin(cst * y);
								Qold[nd * 6 + 3] = A * (cos(cst * x) - cos(cst * z));
								Qold[nd * 6 + 4] = A * sin(cst * x);
								Qold[nd * 6 + 5] = A * (cos(cst * y) - cos(cst * x));

								nd++;
							}
							else if(init_bulktype[l] == 2 || init_bulktype[l] == 3){
								nd++;
							}
							l++;
						}
					}
				}
			}
			else if (seed == -1446){

				for(int k = 0; k < Nz; k++){
					for(int j = 0; j < Ny; j++){
						for(int i = 0; i < Nx; i++){
							
							if(init_bulktype[l] == 1){
								
								xi = (i - rx) * cst * isq2;
								yi = (j - ry) * cst * isq2;
								zi = (k - rz) * cst * isq2;
								
								x = xi;
								y = cos(theta) * yi + sin(theta) * zi;
								z = -sin(theta) * yi + cos(theta) * zi;
								
								Qold[nd * 6 + 0] = A * (cos(cst * z) - cos(cst * y));
								Qold[nd * 6 + 1] = A * sin(cst * z);
								Qold[nd * 6 + 2] = A * sin(cst * y);
								Qold[nd * 6 + 3] = A * (cos(cst * x) - cos(cst * z));
								Qold[nd * 6 + 4] = A * sin(cst * x);
								Qold[nd * 6 + 5] = A * (cos(cst * y) - cos(cst * x));

								Qold[nd * 6 + 0] = A * (- sin(y) * cos(x) - sin(x) * cos(z) + 2 * sin(z) * cos(y));
								Qold[nd * 6 + 3] = A * (- sin(z) * cos(y) - sin(y) * cos(x) + 2 * sin(x) * cos(z));
								Qold[nd * 6 + 5] = A * (- sin(x) * cos(z) - sin(z) * cos(y) + 2 * sin(y) * cos(x));
								Qold[nd * 6 + 1] = A * (- sq2 * sin(x) * sin(z) - sq2 * cos(y) * cos(z) + sin(x) * cos(y));
								Qold[nd * 6 + 2] = A * (- sq2 * sin(z) * sin(y) - sq2 * cos(x) * cos(y) + sin(z) * cos(x));
								Qold[nd * 6 + 4] = A * (- sq2 * sin(y) * sin(x) - sq2 * cos(z) * cos(x) + sin(y) * cos(z));

								nd++;
							}
							else if(init_bulktype[l] == 2 || init_bulktype[l] == 3){
								nd++;
							}
							l++;
						}
					}
				}
			}
		}
	}
	//BPI (seed = 4) and BPII (seed = 5)
	else if(seed == 4 || seed == 5){

		double A = 0.2;
		double cst;
		double isq2 = 1.0 / sqrt(2.);
		double sq2 = sqrt(2.);
					
		cst = 2. * qch * redshift;
			
		l = 0;
		nd = 0;

		if(seed == 4){
			for(int k = 0; k < Nz; k++){
				for (int j = 0; j < Ny; j++){
					for (int i = 0; i < Nx; i++){
						
						if(drop[l] || boundary[l] || nboundary[l]){
							if(geo == -2){
								x = (i - rx) * cst * isq2;
								y = (j - 2) * cst * isq2;
								z = (k) * cst * isq2;
							}
							else if(geo == -3){
								x = (i - rx) * cst * isq2;
								y = (j - 2) * cst * isq2;
								z = (k - rz) * cst * isq2;
							}
							else{
								x = (double)(i - rx) * cst * isq2;
								y = (double)(j - ry) * cst * isq2;
								z = (double)(k - rz) * cst * isq2;
							}
										
							Qold[nd * 6 + 0] = A * (- sin(y) * cos(x) - sin(x) * cos(z) + 2 * sin(z) * cos(y));
							Qold[nd * 6 + 3] = A * (- sin(z) * cos(y) - sin(y) * cos(x) + 2 * sin(x) * cos(z));
							Qold[nd * 6 + 5] = A * (- sin(x) * cos(z) - sin(z) * cos(y) + 2 * sin(y) * cos(x));
							Qold[nd * 6 + 1] = A * (- sq2 * sin(x) * sin(z) - sq2 * cos(y) * cos(z) + sin(x) * cos(y));
							Qold[nd * 6 + 2] = A * (- sq2 * sin(z) * sin(y) - sq2 * cos(x) * cos(y) + sin(z) * cos(x));
							Qold[nd * 6 + 4] = A * (- sq2 * sin(y) * sin(x) - sq2 * cos(z) * cos(x) + sin(y) * cos(z));

							nd ++;
						}
						l ++;

					}
				}
			}
		}
		else{

			//time_t t;
    		//srand((unsigned) time(&t));
			srand(rand_seed);
			double dir_temp[3] = { 0. };

			for(int k = 0; k < Nz; k++){
				for (int j = 0; j < Ny; j++){
					for (int i = 0; i < Nx; i++){
						
						if(drop[l] || boundary[l] || nboundary[l]){
							
							if(geo == -2){

								x = i - rx;
								y = j - 2;
								z = k;

							}
							else if(geo == -3){

								x = i - rx;
								y = j - 2;
								z = k - rz;

							}
							else{

								x = (double)(i - rx) * dx;
								y = (double)(j - ry) * dy;
								z = (double)(k - rz) * dz;

							}

							if(interface != 0 && geo == 10){
								if(init_bulktype[l] == 3){

									dir_temp[0] = (rand() % (pRx + interface) + 1);
						        	dir_temp[1] = (rand() % (pRy + interface) + 1);
						        	dir_temp[2] = (rand() % (pRz + interface) + 1);
                               		norm_v(dir_temp);

									Qold[nd * 6 + 0] = dir2ten(dir_temp, 0, S2);
									Qold[nd * 6 + 1] = dir2ten(dir_temp, 1, S2);
									Qold[nd * 6 + 2] = dir2ten(dir_temp, 2, S2);
									Qold[nd * 6 + 3] = dir2ten(dir_temp, 3, S2);
									Qold[nd * 6 + 4] = dir2ten(dir_temp, 4, S2);
									Qold[nd * 6 + 5] = dir2ten(dir_temp, 5, S2);
								}
								else{
									Qold[nd * 6 + 0] = A * (cos(cst * z) - cos(cst * y));
									Qold[nd * 6 + 1] = A * sin(cst * z);
									Qold[nd * 6 + 2] = A * sin(cst * y);
									Qold[nd * 6 + 3] = A * (cos(cst * x) - cos(cst * z));
									Qold[nd * 6 + 4] = A * sin(cst * x);
									Qold[nd * 6 + 5] = A * (cos(cst * y) - cos(cst * x));
								}
							}
							else{
								Qold[nd * 6 + 0] = A * (cos(cst * z) - cos(cst * y));
								Qold[nd * 6 + 1] = A * sin(cst * z);
								Qold[nd * 6 + 2] = A * sin(cst * y);
								Qold[nd * 6 + 3] = A * (cos(cst * x) - cos(cst * z));
								Qold[nd * 6 + 4] = A * sin(cst * x);
								Qold[nd * 6 + 5] = A * (cos(cst * y) - cos(cst * x));
								
							}
							nd ++;
						}
						l ++;
					}
				}
			}
		}
	}              

	
	//seed=6; [110] BPI; 7: [110] BPII; 8: [111] BPI; 9: [111] BPII
	else if(seed == 6 || seed == 7 || seed==8 || seed==9){
                double A = 0.2;
                double cst,phi;
                double theta = 45 / 180.0 * M_PI;
                double xj,yj,zj;
                //double theta = 90.0 / 180.0 * M_PI;
                double isq2 = 1.0 / sqrt(2);
                double sq2 = sqrt(2);

               
                cst = 2 * qch * redshift;
                
               

                l = 0;
				nd = 0;

                for(k = 0; k < Nz; k++){
                        for (j = 0; j < Ny; j++){
                                for (i = 0; i < Nx; i++){
                                        if(drop[l] || boundary[l] || nboundary[l]){
                                                if(seed == 6){

													if(geo == -2){
														xi = (i - rx) * cst * isq2;
                                                    	yi = (j - 2) * cst * isq2;
                                                   	 	zi = k * cst * isq2;	
													}

													else if(geo == -3){

														xi = (i - rx) * cst * isq2;
                                                    	yi = (j - 2) * cst * isq2;
                                                   	 	zi = (k - ry) * cst * isq2;

													}

													else{

														xi = (i - rx) * cst * isq2;
                                                    	yi = (j - ry) * cst * isq2;
                                                   	 	zi = (k - rz) * cst * isq2;

													}
													
													x = xi;
													y = cos(theta) * yi + sin(theta) * zi;
													z = -sin(theta) * yi + cos(theta) * zi;

                                                    Qold[nd * 6 + 0] = A * (- sin(y) * cos(x) - sin(x) * cos(z) + 2 * sin(z) * cos(y));
                                                    Qold[nd * 6 + 3] = A * (- sin(z) * cos(y) - sin(y) * cos(x) + 2 * sin(x) * cos(z));
                                                    Qold[nd * 6 + 5] = A * (- sin(x) * cos(z) - sin(z) * cos(y) + 2 * sin(y) * cos(x));
                                                    Qold[nd * 6 + 1] = A * (- sq2 * sin(x) * sin(z) - sq2 * cos(y) * cos(z) + sin(x) * cos(y));
                                                    Qold[nd * 6 + 2] = A * (- sq2 * sin(z) * sin(y) - sq2 * cos(x) * cos(y) + sin(z) * cos(x));
                                                    Qold[nd * 6 + 4] = A * (- sq2 * sin(y) * sin(x) - sq2 * cos(z) * cos(x) + sin(y) * cos(z));
                                                }
                                                
												else if(seed == 7){
                                                        xi = i - rx;
                                                        yi = j - ry;
                                                        zi = k - rz;
														x = xi;
														y = cos(theta) * yi + sin(theta) * zi;
														z = -sin(theta) * yi + cos(theta) * zi;
                                                        Qold[nd * 6 + 0] = A * (cos(cst * z) - cos(cst * y));
                                                        Qold[nd * 6 + 1] = A * sin(cst * z);
                                                        Qold[nd * 6 + 2] = A * sin(cst * y);
                                                        Qold[nd * 6 + 3] = A * (cos(cst * x) - cos(cst * z));
                                                        Qold[nd * 6 + 4] = A * sin(cst * x);
                                                        Qold[nd * 6 + 5] = A * (cos(cst * y) - cos(cst * x));
                                                }
						

						
                                                else if(seed == 8){    //8 and 9, (111) planes oriented

                                                        xi = (i - rx) * cst * isq2;
                                                        yi = (j - ry) * cst * isq2;
                                                        zi = (k - rz) * cst * isq2;
                                                       
														theta=atan(1.0/sqrt(2.0));
														//BPI_(211)
                                                        //Rotation around vector (-1,1,0)
                                                       	x=xi*0.5*(1.0+cos(theta))-0.5*yi*(1.0-cos(theta))+zi*sin(theta)/sqrt(2.0); 
                                                      	y=-0.5*xi*(1.0-cos(theta))+yi*0.5*(1.0+cos(theta))+zi*sin(theta)/sqrt(2.0);
                                                       	z=-xi*sin(theta)/sqrt(2.0)-yi*sin(theta)/sqrt(2.0)+zi*cos(theta);
							
                                                       	theta= 1.0*M_PI/12.0;
                                                       	xi=x;
                                                       	yi=y;
                                                       	zi=z;
                                                       	//Rotation around the vector (1,1,1) here it should be implemented rotation around (2,1,1)
                                                       	x = xi*1.0/3.0*(2.0*cos(theta)+1.0) + yi*((1.0-cos(theta))/3.0-sin(theta)/sqrt(3.0)) + zi*((1.0-cos(theta))/3.0+sin(theta)/sqrt(3.0)) ;
                                                       	y = xi*((1.0-cos(theta))/3.0+sin(theta)/sqrt(3.0)) + yi*1.0/3.0*(2.0*cos(theta)+1.0) + zi*((1.0-cos(theta))/3.0-sin(theta)/sqrt(3.0)) ;
                                                       	z = xi*((1.0-cos(theta))/3.0-sin(theta)/sqrt(3.0)) + yi*((1.0-cos(theta))/3.0+sin(theta)/sqrt(3.0)) + zi*1.0/3.0*(2.0*cos(theta)+1.0) ;
                          
                                                           
                                                        Qold[nd * 6 + 0] = A * (- sin(y) * cos(x) - sin(x) * cos(z) + 2 * sin(z) * cos(y));
                                                        Qold[nd * 6 + 3] = A * (- sin(z) * cos(y) - sin(y) * cos(x) + 2 * sin(x) * cos(z));
                                                        Qold[nd * 6 + 5] = A * (- sin(x) * cos(z) - sin(z) * cos(y) + 2 * sin(y) * cos(x));
                                                        Qold[nd * 6 + 1] = A * (- sq2 * sin(x) * sin(z) - sq2 * cos(y) * cos(z) + sin(x) * cos(y));
                                                        Qold[nd * 6 + 2] = A * (- sq2 * sin(z) * sin(y) - sq2 * cos(x) * cos(y) + sin(z) * cos(x));
                                                        Qold[nd * 6 + 4] = A * (- sq2 * sin(y) * sin(x) - sq2 * cos(z) * cos(x) + sin(y) * cos(z));
                                                }
                                                else if(seed == 9){

													if(geo == -2){
														xi = i - rx;
                                                        yi = j - 2;
                                                        zi = k;
													}
													
													else if(geo == -3){
														xi = i - rx;
                                                        yi = j - 2;
                                                        zi = k - rz;
													}
													else{
														xi = i - rx;
                                                        yi = j - ry;
                                                        zi = k - rz;
													}
                                                                     
                                                    theta=atan(sqrt(2.0));

													//BPII_(111)
                                                    //Rotation around vector (-1,1,0)
                                                    x=xi*0.5*(1.0+cos(theta))-0.5*yi*(1.0-cos(theta))+zi*sin(theta)/sqrt(2.0); 
                                                    y=-0.5*xi*(1.0-cos(theta))+yi*0.5*(1.0+cos(theta))+zi*sin(theta)/sqrt(2.0);
                                                    z=-xi*sin(theta)/sqrt(2.0)-yi*sin(theta)/sqrt(2.0)+zi*cos(theta);
							
                                                    theta= 1.0*M_PI/12.0;
                                                    xi=x;
                                                    yi=y;
                                                    zi=z;

                                                    //Rotation around the vector (1,1,1) 
                                                    x = xi*1.0/3.0*(2.0*cos(theta)+1.0) + yi*((1.0-cos(theta))/3.0-sin(theta)/sqrt(3.0)) + zi*((1.0-cos(theta))/3.0+sin(theta)/sqrt(3.0)) ;
                                                    y = xi*((1.0-cos(theta))/3.0+sin(theta)/sqrt(3.0)) + yi*1.0/3.0*(2.0*cos(theta)+1.0) + zi*((1.0-cos(theta))/3.0-sin(theta)/sqrt(3.0)) ;
                                                    z = xi*((1.0-cos(theta))/3.0-sin(theta)/sqrt(3.0)) + yi*((1.0-cos(theta))/3.0+sin(theta)/sqrt(3.0)) + zi*1.0/3.0*(2.0*cos(theta)+1.0) ;
                                                        
													Qold[nd * 6 + 0] = A * (cos(cst * z) - cos(cst * y));
                                                    Qold[nd * 6 + 1] = A * sin(cst * z);
                                                    Qold[nd * 6 + 2] = A * sin(cst * y);
                                                    Qold[nd * 6 + 3] = A * (cos(cst * x) - cos(cst * z));
                                                    Qold[nd * 6 + 4] = A * sin(cst * x);
                                                    Qold[nd * 6 + 5] = A * (cos(cst * y) - cos(cst * x));

                                                }
										nd ++;
                                        }
                                        l ++;

                                }
                        }
                }
        }

	//twist bipolar particle
	else if(seed == 12){
		l = 0;
		nd = 0;
		if(!norm_v(init_dir))   return false;
		for(k = 0; k < Nz; k++){
			for (j = 0; j < Ny; j++){
				for (i = 0; i < Nx; i++){
					if(drop[l] || boundary[l] || nboundary[l]){
						dir[0] = init_dir[0];			
						dir[1] = init_dir[1];			
						dir[2] = init_dir[2];			
						for(n = 0; n < Np; n++){
							x = i - pos[n][0];
							y = j - pos[n][1];
							z = k - pos[n][2];
							if(x > 0.5 * Nx){
								x -= Nx;
							}
							else if(x < -0.5 * Nx){
								x += Nx;
							}
							if(y > 0.5 * Ny){
								y -= Ny;
							}
							else if(y < -0.5 * Ny){
								y += Ny;
							}
							disxy = sqrt(x*x+y*y);
							dis = sqrt(x*x+y*y + z*z);
							if(fabs(z) > Rp && disxy < Rp && disxy > 0 && fabs(z) < 2 * Rp){
								if(n == 0){
									dir[0] += - 2 * z / fabs(z) *  Rp * Rp / dis/ dis * y / disxy;			
									dir[1] += 2 * z / fabs(z) * Rp * Rp / dis/ dis * x / disxy;			
								}
								else{
									dir[0] += 2 * z / fabs(z) * Rp * Rp / dis/ dis * y / disxy;			
									dir[1] += -2 * z / fabs(z) * Rp * Rp / dis/ dis * x / disxy;			
								}
							}
						}
						if(!norm_v(dir))   return false;
						for(n = 0; n < 6; n ++){
							Qold[nd * 6 + n] = dir2ten(dir, n, S);
						}
						nd ++;
					}
					l ++;
				}
			}
		}
        }
	
	//twist bipolar particle
	else if(seed == 22){
		l = 0;
		nd = 0;
		if(!norm_v(init_dir))   return false;
		for(k = 0; k < Nz; k++){
			for (j = 0; j < Ny; j++){
				for (i = 0; i < Nx; i++){
					if(drop[l] || boundary[l] || nboundary[l]){
						dir[0] = init_dir[0];			
						dir[1] = init_dir[1];			
						dir[2] = init_dir[2];			
						for(n = 0; n < Np; n++){
							x = i - pos[n][0];
							y = j - pos[n][1];
							z = k - pos[n][2];
							if(x > 0.5 * Nx){
								x -= Nx;
							}
							else if(x < -0.5 * Nx){
								x += Nx;
							}
							if(y > 0.5 * Ny){
								y -= Ny;
							}
							else if(y < -0.5 * Ny){
								y += Ny;
							}
							disxy = sqrt(x*x+y*y);
							dis = sqrt(x*x+y*y + z*z);
							if(fabs(z) > Rp && disxy < Rp && disxy > 0 && fabs(z) < 2 * Rp){
								if(n == 0){
									dir[0] += - 2 * Rp * Rp / dis/ dis * y / disxy;			
									dir[1] += 2 * Rp * Rp / dis/ dis * x / disxy;			
								}
								else{
									dir[0] += 2 * Rp * Rp / dis/ dis * y / disxy;			
									dir[1] += -2 * Rp * Rp / dis/ dis * x / disxy;			
								}
							}
						}
						if(!norm_v(dir))   return false;
						for(n = 0; n < 6; n ++){
							Qold[nd * 6 + n] = dir2ten(dir, n, S);
						}
						nd ++;
					}
					l ++;
				}
			}
		}
        }
	
	//twist bipolar particle
	else if(seed == 13){
		l = 0;
		nd = 0;
		if(!norm_v(init_dir))   return false;
		for(k = 0; k < Nz; k++){
			for (j = 0; j < Ny; j++){
				for (i = 0; i < Nx; i++){
					if(drop[l] || boundary[l] || nboundary[l]){
						dir[0] = init_dir[0];			
						dir[1] = init_dir[1];			
						dir[2] = init_dir[2];			
						for(n = 0; n < Np; n++){
							x = i - pos[n][0];
							y = j - pos[n][1];
							z = k - pos[n][2];
							if(x > 0.5 * Nx){
								x -= Nx;
							}
							else if(x < -0.5 * Nx){
								x += Nx;
							}
							if(y > 0.5 * Ny){
								y -= Ny;
							}
							else if(y < -0.5 * Ny){
								y += Ny;
							}
							disxy = sqrt(x*x+y*y);
							dis = sqrt(x*x+y*y + z*z);
							if(fabs(z) > Rp && disxy < Rp && disxy > 0 && fabs(z) < 2 * Rp){
								dir[0] += - 2 * z / fabs(z) * Rp * Rp / dis/ dis * y / disxy;			
								dir[1] += 2 * z / fabs(z) * Rp * Rp / dis/ dis * x / disxy;			
							}
						}
						if(!norm_v(dir))   return false;
						for(n = 0; n < 6; n ++){
							Qold[nd * 6 + n] = dir2ten(dir, n, S);
						}
						nd ++;
					}
					l ++;
				}
			}
		}
        }
	
	//twist bipolar particle
	else if(seed == 23){
		l = 0;
		nd = 0;
		if(!norm_v(init_dir))   return false;
		for(k = 0; k < Nz; k++){
			for (j = 0; j < Ny; j++){
				for (i = 0; i < Nx; i++){
					if(drop[l] || boundary[l] || nboundary[l]){
						dir[0] = init_dir[0];			
						dir[1] = init_dir[1];			
						dir[2] = init_dir[2];			
						for(n = 0; n < Np; n++){
							x = i - pos[n][0];
							y = j - pos[n][1];
							z = k - pos[n][2];
							if(x > 0.5 * Nx){
								x -= Nx;
							}
							else if(x < -0.5 * Nx){
								x += Nx;
							}
							if(y > 0.5 * Ny){
								y -= Ny;
							}
							else if(y < -0.5 * Ny){
								y += Ny;
							}
							disxy = sqrt(x*x+y*y);
							dis = sqrt(x*x+y*y + z*z);
							if(fabs(z) > Rp && disxy < Rp && disxy > 0 && fabs(z) < 2 * Rp){
									dir[0] += - 2 * Rp * Rp / dis/ dis * y / disxy;			
									dir[1] += 2 * Rp * Rp / dis/ dis * x / disxy;			
							}
						}
						if(!norm_v(dir))   return false;
						for(n = 0; n < 6; n ++){
							Qold[nd * 6 + n] = dir2ten(dir, n, S);
						}
						nd ++;
					}
					l ++;
				}
			}
		}
        }
	
	//twist bipolar particle
	else if(seed == 14){
		l = 0;
		nd = 0;
		if(!norm_v(init_dir))   return false;
		for(k = 0; k < Nz; k++){
			for (j = 0; j < Ny; j++){
				for (i = 0; i < Nx; i++){
					if(drop[l] || boundary[l] || nboundary[l]){
						dir[0] = init_dir[0];			
						dir[1] = init_dir[1];			
						dir[2] = init_dir[2];			
						for(n = 0; n < Np; n++){
							x = i - pos[n][0];
							y = j - pos[n][1];
							z = k - pos[n][2];
							if(x > 0.5 * Nx){
								x -= Nx;
							}
							else if(x < -0.5 * Nx){
								x += Nx;
							}
							if(y > 0.5 * Ny){
								y -= Ny;
							}
							else if(y < -0.5 * Ny){
								y += Ny;
							}
							dis = sqrt(x*x+y*y+z*z);
							disxy = sqrt(x*x+y*y);
							if(fabs(z) > Rp && disxy > 0){
								dir[0] += 2 * z / fabs(z) * Rp * Rp / dis/ dis * y / disxy;			
								dir[1] += - 2 * z / fabs(z) * Rp * Rp / dis/ dis * x / disxy;			
							}
						}
						if(!norm_v(dir))   return false;
						for(n = 0; n < 6; n ++){
							Qold[nd * 6 + n] = dir2ten(dir, n, S);
						}
						nd ++;
					}
					l ++;
				}
			}
		}
        }
	
		//radial droplet and axial cylinders.
		//Se han agregado configuraciones iniciales para medios cilindros y medias gotas. Para el medio cilindro y cilindro completo adopta una
		//geometría axial mientras que para la gota y la media gota, adquiere una geometría inicial radial.

		//Axial and Radial configurations.
    else if(seed == 11){
        	l = 0;
			nd = 0;
			double mod = 0;

			if(geo == 2){

				for(k = 0; k < Nz; k++){
            		for (j = 0; j < Ny; j++){
                		for (i = 0; i < Nx; i++){
                   
				   			if(drop[l] || boundary[l] || nboundary[l]){
                        		x = (i - rx) * dx;
                        		y = (j - ry) * dy;//j * dy;

								//// ***** RECORDAR MODIFICAR LA CONDICION DE Z *********////

                        		z = 5 * dz; //10; //k * dz;    //Modifiqué z = (k - rz) * dz por k * dz

								dir[0] = -x;
								dir[1] = -y;
								dir[2] = -z; //-z;
						
								mod = sqrt(dir[0] * dir[0] + dir[1] * dir[1] + dir[2] * dir[2]); 

								if (mod == 0){
									dir[0] = 0;
                            		dir[1] = 0;
                            		dir[2] = -1; // -10;
									//printf("i: %d, j: %d k: %d \n", i, j, k);
								}

								else{
									dir[0] = dir[0] / mod;
									dir[1] = dir[1] / mod;
									dir[2] = dir[2] / mod;						
								}												

                        		for (n = 0; n < 6; n++) {
                    	    		Qold[nd * 6 + n] = dir2ten(dir, n, S);
                        		}
								nd ++;
                   			}

							l ++;

                		}
            		}
        		}

			}

			else if(geo == -2){
				for(k = 0; k < Nz; k++){
            		for (j = 0; j < Ny; j++){
                		for (i = 0; i < Nx; i++){
                   
				   			if(drop[l] || boundary[l] || nboundary[l]){
                        		

								if(j > 2 && boundary[l] == true){

									x = (i - rx) * dx;
                        			y = (j - 2) * dy;//j * dy;
                        			z = 0;

									dir[0] = -x;
									dir[1] = -y;
									dir[2] = -z; 

									mod = sqrt(dir[0] * dir[0] + dir[1] * dir[1] + dir[2] * dir[2]); 									

									dir[0] = dir[0] / mod;
                            		dir[1] = dir[1] / mod;
                            		dir[2] = dir[2] / mod;

								}

								else{

									x = (i - rx) * dx;
                        			y = (j - 2) * dy;//j * dy;
                        			z = (5) * dz; //10; //k * dz;    //Modifiqué z = (k - rz) * dz por k * dz

									dir[0] = -x;
									dir[1] = -y;
									dir[2] = -z; //-z;
						
									mod = sqrt(dir[0] * dir[0] + dir[1] * dir[1] + dir[2] * dir[2]); 

									if (mod == 0){
										dir[0] = 0;
                            			dir[1] = 0;
                            			dir[2] = -1; // -10;
										//printf("i: %d, j: %d k: %d \n", i, j, k);
									}

									else{
										dir[0] = dir[0] / mod;
										dir[1] = dir[1] / mod;
										dir[2] = dir[2] / mod;						
									}

								}																			

                        		for (n = 0; n < 6; n++) {
                    	    		Qold[nd * 6 + n] = dir2ten(dir, n, S);
                        		}
								nd ++;
                   			}

							l ++;

                		}
            		}
        		}
			}

			else if(geo == -22){
				for(k = 0; k < Nz; k++){
            		for (j = 0; j < Ny; j++){
                		for (i = 0; i < Nx; i++){
                   
				   			if(drop[l] || boundary[l] || nboundary[l]){
                        		

								if(j > 2 && boundary[l] == true){

									x = (i - 2) * dx;
                        			y = (j - 2) * dy;//j * dy;
                        			z = 0;

									dir[0] = -x;
									dir[1] = -y;
									dir[2] = -z; 

									mod = sqrt(dir[0] * dir[0] + dir[1] * dir[1] + dir[2] * dir[2]); 									

									dir[0] = dir[0] / mod;
                            		dir[1] = dir[1] / mod;
                            		dir[2] = dir[2] / mod;

								}

								else{

									x = (i - 2) * dx;
                        			y = (j - 2) * dy;//j * dy;
                        			z = (5) * dz; //10; //k * dz;    //Modifiqué z = (k - rz) * dz por k * dz

									dir[0] = -x;
									dir[1] = -y;
									dir[2] = -z; //-z;
						
									mod = sqrt(dir[0] * dir[0] + dir[1] * dir[1] + dir[2] * dir[2]); 

									if (mod == 0){
										dir[0] = 0;
                            			dir[1] = 0;
                            			dir[2] = -1; // -10;
										//printf("i: %d, j: %d k: %d \n", i, j, k);
									}

									else{
										dir[0] = dir[0] / mod;
										dir[1] = dir[1] / mod;
										dir[2] = dir[2] / mod;						
									}

								}																			

                        		for (n = 0; n < 6; n++) {
                    	    		Qold[nd * 6 + n] = dir2ten(dir, n, S);
                        		}
								nd ++;
                   			}

							l ++;

                		}
            		}
        		}
			}

			else if(geo == 3){

				for(k = 0; k < Nz; k++){
            		for (j = 0; j < Ny; j++){
                		for (i = 0; i < Nx; i++){
                   
				   			if(drop[l] || boundary[l] || nboundary[l]){
                        		x = (i - rx) * dx;
                        		y = (j - ry) * dy;//j * dy;
                        		z = (k - rz) * dz; //10; //k * dz;    //Modifiqué z = (k - rz) * dz por k * dz

								dir[0] = -x;
								dir[1] = -y;
								dir[2] = -z; //-z;
						
								mod = sqrt(dir[0] * dir[0] + dir[1] * dir[1] + dir[2] * dir[2]); 

								if (mod == 0){
									dir[0] = 0;
                            		dir[1] = 0;
                            		dir[2] = -1; // -10;
									//printf("i: %d, j: %d k: %d \n", i, j, k);
								}

								else{
									dir[0] = dir[0] / mod;
									dir[1] = dir[1] / mod;
									dir[2] = dir[2] / mod;						
								}												

                        		for (n = 0; n < 6; n++) {
                    	    		Qold[nd * 6 + n] = dir2ten(dir, n, S);
                        		}
								nd ++;
                   			}

							l ++;

                		}
            		}
        		}

			}

			//radial para media esfera
			else if(geo == -3){

				for(k = 0; k < Nz; k++){
            		for (j = 0; j < Ny; j++){
                		for (i = 0; i < Nx; i++){
                   
				   			if(drop[l] || boundary[l] || nboundary[l]){
                        		x = (i - rx) * dx;
                        		y = (j - 2) * dy;//j * dy;
                        		z = (k - rz) * dz; //10; //k * dz;    //Modifiqué z = (k - rz) * dz por k * dz

								dir[0] = -x;
								dir[1] = -y;
								dir[2] = -z; //-z;
						
								mod = sqrt(dir[0] * dir[0] + dir[1] * dir[1] + dir[2] * dir[2]); 

								if (mod == 0){
									dir[0] = 0;
                            		dir[1] = 0;
                            		dir[2] = -1; // -10;
									//printf("i: %d, j: %d k: %d \n", i, j, k);
								}

								else{
									dir[0] = dir[0] / mod;
									dir[1] = dir[1] / mod;
									dir[2] = dir[2] / mod;						
								}												

                        		for (n = 0; n < 6; n++) {
                    	    		Qold[nd * 6 + n] = dir2ten(dir, n, S);
                        		}
								nd ++;
                   			}

							l ++;

                		}
            		}
        		}

			}

    		//radial para 1/4 de esfera
			else if(geo == -33){	
				for(k = 0; k < Nz; k++){
            		for (j = 0; j < Ny; j++){
                		for (i = 0; i < Nx; i++){
                   
				   			if(drop[l] || boundary[l] || nboundary[l]){
                        		x = (i - 2) * dx;	//aquí es por el espaciamiento en x debido a que es 1/4 de esfera
                        		y = (j - 2) * dy; 	
                        		z = (k - rz) * dz;

								dir[0] = -x;
								dir[1] = -y;
								dir[2] = -z; //-z;
						
								mod = sqrt(dir[0] * dir[0] + dir[1] * dir[1] + dir[2] * dir[2]); 

								if (mod == 0){
									dir[0] = 0;
                            		dir[1] = 0;
                            		dir[2] = -1; // -10;
									//printf("i: %d, j: %d k: %d \n", i, j, k);
								}

								else{
									dir[0] = dir[0] / mod;
									dir[1] = dir[1] / mod;
									dir[2] = dir[2] / mod;						
								}												

                        		for (n = 0; n < 6; n++) {
                    	    		Qold[nd * 6 + n] = dir2ten(dir, n, S);
                        		}
								nd ++;
                   			}

							l ++;

                		}
            		}
        		}
			}

			//Radial para elipsoide con dos U.
			else if(geo == 4){

				for(int k = 0; k < Nz; k++){
            		for (int j = 0; j < Ny; j++){
                		for (int i = 0; i < Nx; i++){
                   
				   			if(drop[l] || boundary[l] || nboundary[l]){
                        		x = (i - rx) * dx;	//aquí es por el espaciamiento en x debido a que es 1/4 de esfera
                        		y = (j - ry) * dy; 	
                        		z = (k - rz) * dz;

								dir[0] = -x;
								dir[1] = -y;
								dir[2] = -z; //-z;
						
								mod = sqrt(dir[0] * dir[0] + dir[1] * dir[1] + dir[2] * dir[2]); 

								if (mod == 0){
									dir[0] = 0;
                            		dir[1] = 0;
                            		dir[2] = -1; // -10;
									//printf("i: %d, j: %d k: %d \n", i, j, k);
								}

								else{
									dir[0] = dir[0] / mod;
									dir[1] = dir[1] / mod;
									dir[2] = dir[2] / mod;						
								}

								if(DoubleU){
									if(bulktype[l] == 1){

										for (n = 0; n < 6; n++) {
											Qold[nd * 6 + n] = dir2ten(dir, n, S);
										}

									}
									else if(bulktype[l] == 2){

										for (n = 0; n < 6; n++) {
											Qold[nd * 6 + n] = dir2ten(dir, n, S2);
										}

									}
								}
								else{

									for (n = 0; n < 6; n++) {
										Qold[nd * 6 + n] = dir2ten(dir, n, S);
									}

								}
								
								nd ++;
                   			}

							l ++;

                		}
            		}
        		}

			}
		}

	//escaped cylinder
	else if(seed == 21){
		l = 0;
		nd = 0;
		double Rx = rx - 2;
		double Ry = ry - 2;
		for(k = 0; k < Nz; k++){
			for (j = 0; j < Ny; j++){
				for (i = 0; i < Nx; i++){
					if(drop[l] ){//|| boundary[l]){ //|| nboundary[l]){
						x = i - rx;
						y = j - ry;
						dis = sqrt(x*x+(y+Ry)*(y+Ry));
						dir[0] = x / Rx;
						dir[1] = y / Ry;
						dir[2] = sqrt((Ry * Rx - dis * dis) / (Ry*Rx));
					//	dis=sqrt(dir[0]*dir[0]+dir[1]*dir[1]+ dir[2]* dir[2]);
					//	 dir[0] =  dir[0] /dis;
					//	 dir[1] =  dir[1] /dis;
					//	 dir[2] =  dir[2] /dis;
						if(!norm_v(dir))   return false;
						for(n = 0; n < 6; n ++){
							Qold[nd * 6 + n] = dir2ten(dir, n, S);
						}
						nd ++;
					}
					l ++;
				}
			}
		}
        }
	//helical along z
    else if(seed == 87){
                l = 0;
		nd = 0;
                for(k = 0; k < Nz; k++){
                        for (j = 0; j < Ny; j++){
                                for (i = 0; i < Nx; i++){
                                        if(drop[l] || boundary[l] || nboundary[l]){
						dir[0] = cos(qch * (k - rz));
						dir[1] = sin(qch * (k - rz));
						for(n = 0; n < 6; n ++){
							Qold[nd * 6 + n] = dir2ten(dir, n, S);
						}
						nd ++;
                                        }
                                        l ++;
                                }
                        }
                }
    }
	//helical along x
    else if(seed == 88){
                l = 0;
		nd = 0;
                for(k = 0; k < Nz; k++){
                        for (j = 0; j < Ny; j++){
                                for (i = 0; i < Nx; i++){
                                        if(drop[l] || boundary[l] || nboundary[l]){
						dir[1] = cos(qch * (i - rx));
						dir[2] = sin(qch * (i - rx));
						for(n = 0; n < 6; n ++){
							Qold[nd * 6 + n] = dir2ten(dir, n, S);
						}
						nd ++;
                                        }
                                        l ++;
                                }
                        }
                }
        }
	//helical along y
    else if(seed == 89){
        l = 0;
		nd = 0;
			for(k = 0; k < Nz; k++){
				for (j = 0; j < Ny; j++){
						for (i = 0; i < Nx; i++){
								if(drop[l] || boundary[l] || nboundary[l]){
				dir[0] = cos(qch * (j - ry));
				dir[2] = sin(qch * (j - ry));
				for(n = 0; n < 6; n ++){
					Qold[nd * 6 + n] = dir2ten(dir, n, S);
				}
				nd ++;
								}
								l ++;
						}
				}
			}
        }

	
	// ********** 11 + x = condición radial + condición BP *********** //////////////
	// ***** For BPI seed = 114 and 116 with [110]; and for BPII seed = 115 and 119 with [111] ***** //////////
	// ***** seeds 124, 126, 125, 129 and so on are combinations (with a 1 at the beginning) of seed 2 o 3 (cholesteric)
	// ***** and BPI's seeds 4 & 6, and BPII's seeds 5 & 9.
	// *** for seed 141, 142, 143 BPII 100, 110, 111 outter shells with BPIII random inner core //

	else if(seed == 114 || seed == 116 || seed == 115 || seed == 119 ||	
			seed == 124 || seed == 126 || seed == 125 || seed == 129 ||
			seed == 134 || seed == 136 || seed == 135 || seed == 139 || 
			seed == 141 || seed == 142 || seed == 143 || 
			seed == 874 || seed == 875 || seed == 876 || seed == 879)
			{

		double A = 0.2;
		double cst, phi;
		double isq2 = 1.0 / sqrt(2);
		double sq2 = sqrt(2);
        double theta = 45 / 180.0 * M_PI;
        double xj,yj,zj;
		double mod;
		srand(rand_seed);
		double norm = 0.0;
		double dirvec1[3] = {0};

		if( seed == 114 || seed == 116 || seed == 124 || seed == 126 || seed == 134 || seed == 136 || 
			seed == 141 || seed == 142 || seed == 143 ||
			seed == 874 || seed == 876){
			cst = 2 * qch * 0.71;
		}        
		else{
			cst = 2 * qch * 0.86;
		}	

		l = 0;
		nd = 0;

		for(k = 0; k < Nz; k++){
			for (j = 0; j < Ny; j++){
				for (i = 0; i < Nx; i++){

					if(init_bulktype[l] == 1){
						if(drop[l] || boundary[l] || nboundary[l]){
							if(seed == 114 || seed == 116 || seed == 115 || seed == 119){
								
								
								x = (i - rx) * dx;
								y = (j - ry) * dy;//j * dy;
								z = (k - rz) * dz; //10; //k * dz;    //Modifiqué z = (k - rz) * dz por k * dz

								dir[0] = -x;
								dir[1] = -y;
								dir[2] = -z; //-z;
						
								mod = sqrt(dir[0] * dir[0] + dir[1] * dir[1] + dir[2] * dir[2]); 

								if (mod == 0){
									dir[0] = 0;
									dir[1] = 0;
									dir[2] = -1; // -10;
									//printf("i: %d, j: %d k: %d \n", i, j, k);
								}

								else{
									dir[0] = dir[0] / mod;
									dir[1] = dir[1] / mod;
									dir[2] = dir[2] / mod;						
								}												

								for (n = 0; n < 6; n++) {
									Qold[nd * 6 + n] = dir2ten(dir, n, S);
								}
								nd ++;
							
							}
							else if(seed == 141 || seed == 142 || seed == 143){
								
									for(n = 0; n < 3; n ++){
										dir[n] = (double)rand() / (double)RAND_MAX * 2 - 1;
									}

									if(!norm_v(dir)){
										printf("Problems with random initialization\n");
										return false;
									}        	

									for(n = 0; n < 6; n ++){
										Qold[nd * 6 + n] = dir2ten(dir, n, 0.5);
									}            
								 nd ++;
							}
							else if(seed == 124 || seed == 126 || seed == 125 || seed == 129 ||
									seed == 134 || seed == 136 || seed == 135 || seed == 139){
					
								x = (i-rx)*dx;
								y = (j-ry)*dy;
								z = (k-rz)*dz;
								dis = sqrt(x*x+y*y+z*z);

								if(seed == 124 || seed == 126 || seed == 125 || seed == 129){
									omega = dis * qch;
								}   
								
								else{
									omega = atan2(y, x) + dis * qch;
									disxy = sqrt(x * x + y * y);
								}
									
								if(disxy == 0){
									dir[2] = 1;
									dir[0] = dir[1] = 0;
								}

								else{
									costhe = z / dis;
									sinthe = disxy / dis;
									cosphi = x / disxy;
									sinphi = y / disxy;
									dir[0] = cos(omega) * costhe * cosphi - sin(omega) * sinphi;
									dir[1] = cos(omega) * costhe * sinphi + sin(omega) * cosphi;
									dir[2] = - cos(omega) * sinthe;
								}
				
								if(!norm_v(dir)){
									return false;
								}
									
								for (n = 0; n < 6; n++) {
									Qold[nd * 6 + n] = dir2ten(dir, n, S);
								}
						
								nd ++;
							}
							else if(seed == 874 || seed == 875 || seed == 876 || seed == 879){
								
								dir[0] = cos(qch * (k - rz));
								dir[1] = sin(qch * (k - rz));

								for(int n = 0; n < 6; n ++){
									Qold[nd * 6 + n] = dir2ten(dir, n, S);
								}
								nd ++;								
							}
						}
					}

					else if(init_bulktype[l] == 2){
					
						if(drop[l] || boundary[l] || nboundary[l]){
							if(seed == 114 || seed == 124 || seed == 134 || seed == 874){

								x = (i - rx) * cst * isq2;
								y = (j - ry) * cst * isq2;
								z = (k - rz) * cst * isq2;
						
								Qold[nd * 6 + 0] = A * (- sin(y) * cos(x) - sin(x) * cos(z) + 2 * sin(z) * cos(y));
								Qold[nd * 6 + 3] = A * (- sin(z) * cos(y) - sin(y) * cos(x) + 2 * sin(x) * cos(z));
								Qold[nd * 6 + 5] = A * (- sin(x) * cos(z) - sin(z) * cos(y) + 2 * sin(y) * cos(x));
								Qold[nd * 6 + 1] = A * (- sq2 * sin(x) * sin(z) - sq2 * cos(y) * cos(z) + sin(x) * cos(y));
								Qold[nd * 6 + 2] = A * (- sq2 * sin(z) * sin(y) - sq2 * cos(x) * cos(y) + sin(z) * cos(x));
								Qold[nd * 6 + 4] = A * (- sq2 * sin(y) * sin(x) - sq2 * cos(z) * cos(x) + sin(y) * cos(z));
							}

							else if(seed == 115 || seed == 125 || seed == 135 || seed == 141 || seed == 875){
										
								x = i - rx;
								y = j - ry;
								z = k - rz;

								Qold[nd * 6 + 0] = A * (cos(cst * z) - cos(cst * y));
								Qold[nd * 6 + 1] = A * sin(cst * z);
								Qold[nd * 6 + 2] = A * sin(cst * y);
								Qold[nd * 6 + 3] = A * (cos(cst * x) - cos(cst * z));
								Qold[nd * 6 + 4] = A * sin(cst * x);
								Qold[nd * 6 + 5] = A * (cos(cst * y) - cos(cst * x));
							}

							else if(seed == 116  || seed == 126 || seed == 136 || seed == 876){

								xi = (i - rx) * cst * isq2;
								yi = (j - ry) * cst * isq2;
								zi = (k - rz) * cst * isq2;
								
								x = xi;
								y = cos(theta) * yi + sin(theta) * zi;
								z = -sin(theta) * yi + cos(theta) * zi;

								Qold[nd * 6 + 0] = A * (- sin(y) * cos(x) - sin(x) * cos(z) + 2 * sin(z) * cos(y));
								Qold[nd * 6 + 3] = A * (- sin(z) * cos(y) - sin(y) * cos(x) + 2 * sin(x) * cos(z));
								Qold[nd * 6 + 5] = A * (- sin(x) * cos(z) - sin(z) * cos(y) + 2 * sin(y) * cos(x));
								Qold[nd * 6 + 1] = A * (- sq2 * sin(x) * sin(z) - sq2 * cos(y) * cos(z) + sin(x) * cos(y));
								Qold[nd * 6 + 2] = A * (- sq2 * sin(z) * sin(y) - sq2 * cos(x) * cos(y) + sin(z) * cos(x));
								Qold[nd * 6 + 4] = A * (- sq2 * sin(y) * sin(x) - sq2 * cos(z) * cos(x) + sin(y) * cos(z));
							}
													
							else if (seed == 142){
							
								xi = i - rx;
                        		yi = j - ry;
                            	zi = k - rz;

								x = xi;
								y = cos(theta) * yi + sin(theta) * zi;
								z = -sin(theta) * yi + cos(theta) * zi;

								Qold[nd * 6 + 0] = A * (cos(cst * z) - cos(cst * y));
								Qold[nd * 6 + 1] = A * sin(cst * z);
								Qold[nd * 6 + 2] = A * sin(cst * y);
								Qold[nd * 6 + 3] = A * (cos(cst * x) - cos(cst * z));
								Qold[nd * 6 + 4] = A * sin(cst * x);
								Qold[nd * 6 + 5] = A * (cos(cst * y) - cos(cst * x));}

							else if(seed == 119  || seed == 129 || seed == 139 || seed == 143 || seed == 879){

								xi = i - rx;
								yi = j - ry;
								zi = k - rz;
													
								theta=atan(sqrt(2.0));

								//BPII_(111)
								//Rotation around vector (-1,1,0)
								x=xi*0.5*(1.0+cos(theta))-0.5*yi*(1.0-cos(theta))+zi*sin(theta)/sqrt(2.0); 
								y=-0.5*xi*(1.0-cos(theta))+yi*0.5*(1.0+cos(theta))+zi*sin(theta)/sqrt(2.0);
								z=-xi*sin(theta)/sqrt(2.0)-yi*sin(theta)/sqrt(2.0)+zi*cos(theta);

								theta= 1.0*M_PI/12.0;
								xi=x;
								yi=y;
								zi=z;

								//Rotation around the vector (1,1,1) 
								x = xi*1.0/3.0*(2.0*cos(theta)+1.0) + yi*((1.0-cos(theta))/3.0-sin(theta)/sqrt(3.0)) + zi*((1.0-cos(theta))/3.0+sin(theta)/sqrt(3.0)) ;
								y = xi*((1.0-cos(theta))/3.0+sin(theta)/sqrt(3.0)) + yi*1.0/3.0*(2.0*cos(theta)+1.0) + zi*((1.0-cos(theta))/3.0-sin(theta)/sqrt(3.0)) ;
								z = xi*((1.0-cos(theta))/3.0-sin(theta)/sqrt(3.0)) + yi*((1.0-cos(theta))/3.0+sin(theta)/sqrt(3.0)) + zi*1.0/3.0*(2.0*cos(theta)+1.0) ;
									
								Qold[nd * 6 + 0] = A * (cos(cst * z) - cos(cst * y));
								Qold[nd * 6 + 1] = A * sin(cst * z);
								Qold[nd * 6 + 2] = A * sin(cst * y);
								Qold[nd * 6 + 3] = A * (cos(cst * x) - cos(cst * z));
								Qold[nd * 6 + 4] = A * sin(cst * x);
								Qold[nd * 6 + 5] = A * (cos(cst * y) - cos(cst * x));

							}

							nd ++;
						}
					}
					l ++;
				}
			}
		}
	}


	return true;
}
