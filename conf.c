#include "finite.h" 

bool conf(double **pos){
	int l;
	int nb, nd;
	int i, j, k, n;
	double dis, x = 0.0, y = 0.0, z = 0.0, xi, yi, zi;
	int rx = lrint(Nx *0.5);
	int ry = lrint(Ny *0.5);
	int rz = lrint(Nz *0.5);
	double dx = Lx/(Nx-1);
	double dy = Ly/(Ny-1);
	double dz = Lz/(Nz-1);
	double disxy;
	double sinthe, costhe, sinphi, cosphi, omega;
	double Qini[6] = {0};
	double dir[3] = {0};
	int flag = 0;
	bool ideal = false;

	//uniform initial configuration
	if(vseed == 0){

		printf("Uniform configuration.\n");
		nd = 0;
		double dirvec1[3] = {0};
		double dirvec2[3] = {0};
		double norm = 0;

		if((geo == -2 || geo == -3) && ideal == 1){

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

		else{

			if(!norm_v(init_dir)){			
				printf("Problems in initial direction! \n");
				return false;
			}			

			for(n = 0; n < 6; n ++){
				Qini[n] = dir2ten(init_dir, n, S);
			}

			for(nd = 0; nd < droplet; nd ++){
				for(n = 0; n < 6; n ++){
					Qold[nd * 6 + n] = Qini[n];
				}
			}
		}    	

	}

	else if(vseed == 1){
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
	else if(vseed == 10){
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
	else if(vseed == 2 || vseed == 3){
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

						if(vseed == 2){
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
    else if(vseed == -1){
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

	//BPI (vseed = 4) and BPII (vseed = 5)
	else if(vseed == 4 || vseed == 5){

    double A = 0.2;
    double cst;
    double isq2 = 1.0 / sqrt(2);
    double sq2 = sqrt(2);
                
	if(vseed == 4){
        cst = 2 * qch * 0.71;
    }        
	else{
       	cst = 2 * qch * 0.86;
    }

    l = 0;
	nd = 0;

	for(int k = 0; k < Nz; k++){
        for (int j = 0; j < Ny; j++){
            for (int i = 0; i < Nx; i++){
                
				if(drop[l] || boundary[l] || nboundary[l]){
                	if(vseed == 4){

						if(geo == -2){
							x = (i - Nx * 0.5) * cst * isq2;
                        	y = (j - 2) * cst * isq2;
                        	z = (k) * cst * isq2;
						}
						else if(geo == -3){
							x = (i - Nx * 0.5) * cst * isq2;
                        	y = (j - 2) * cst * isq2;
                        	z = (k - Nz * 0.5) * cst * isq2;
						}
						else{
							x = (i - Nx * 0.5) * cst * isq2;
                        	y = (j - Ny * 0.5) * cst * isq2;
                        	z = (k - Nz * 0.5) * cst * isq2;
						}
                                    
						Qold[nd * 6 + 0] = A * (- sin(y) * cos(x) - sin(x) * cos(z) + 2 * sin(z) * cos(y));
                        Qold[nd * 6 + 3] = A * (- sin(z) * cos(y) - sin(y) * cos(x) + 2 * sin(x) * cos(z));
                        Qold[nd * 6 + 5] = A * (- sin(x) * cos(z) - sin(z) * cos(y) + 2 * sin(y) * cos(x));
                        Qold[nd * 6 + 1] = A * (- sq2 * sin(x) * sin(z) - sq2 * cos(y) * cos(z) + sin(x) * cos(y));
                        Qold[nd * 6 + 2] = A * (- sq2 * sin(z) * sin(y) - sq2 * cos(x) * cos(y) + sin(z) * cos(x));
                        Qold[nd * 6 + 4] = A * (- sq2 * sin(y) * sin(x) - sq2 * cos(z) * cos(x) + sin(y) * cos(z));
                    }

                    else if(vseed == 5){
                                
						if(geo == -2){

							x = i - Nx * 0.5;
                    		y = j - 2;
                   			z = k;

						}
						else if(geo == -3){

							x = i - Nx * 0.5;
                    		y = j - 2;
                   			z = k - Nz * 0.5;

						}
						else{

							x = i - Nx * 0.5;
                    		y = j - Ny * 0.5;
                   			z = k - Nz * 0.5;
							   
						}

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
	
	//vseed=6; [110] BPI; 7: [110] BPII; 8: [111] BPI; 9: [111] BPII
	else if(vseed == 6 || vseed == 7 || vseed==8 || vseed==9){
                double A = 0.2;
                double cst,phi;
                double theta = 45 / 180.0 * M_PI;
                double xj,yj,zj;
                //double theta = 90.0 / 180.0 * M_PI;
                double isq2 = 1.0 / sqrt(2);
                double sq2 = sqrt(2);

                if(vseed == 6 || vseed==8){
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
                                        if(drop[l] || boundary[l] || nboundary[l]){
                                                if(vseed == 6){

													if(geo == -2){
														xi = (i - Nx * 0.5) * cst * isq2;
                                                    	yi = (j - 2) * cst * isq2;
                                                   	 	zi = k * cst * isq2;	
													}

													else if(geo == -3){

														xi = (i - Nx * 0.5) * cst * isq2;
                                                    	yi = (j - 2) * cst * isq2;
                                                   	 	zi = (k - Nz * 0.5) * cst * isq2;

													}

													else{

														xi = (i - Nx * 0.5) * cst * isq2;
                                                    	yi = (j - Ny * 0.5) * cst * isq2;
                                                   	 	zi = (k - Nz * 0.5) * cst * isq2;

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
                                                
												else if(vseed == 7){
                                                        xi = i - Nx * 0.5;
                                                        yi = j - Ny * 0.5;
                                                        zi = k - Nz * 0.5;
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
						

						
                                                else if(vseed == 8){    //8 and 9, (111) planes oriented

                                                        xi = (i - Nx * 0.5) * cst * isq2;
                                                        yi = (j - Ny * 0.5) * cst * isq2;
                                                        zi = (k - Nz * 0.5) * cst * isq2;
                                                       
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
                                                else if(vseed == 9){

													if(geo == -2){
														xi = i - Nx * 0.5;
                                                        yi = j - 2;
                                                        zi = k;
													}
													
													else if(geo == -3){
														xi = i - Nx * 0.5;
                                                        yi = j - 2;
                                                        zi = k - Nz * 0.5;
													}
													else{
														xi = i - Nx * 0.5;
                                                        yi = j - Ny * 0.5;
                                                        zi = k - Nz * 0.5;
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
	else if(vseed == 12){
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
	else if(vseed == 22){
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
	else if(vseed == 13){
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
	else if(vseed == 23){
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
	else if(vseed == 14){
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
    else if(vseed == 11){
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
	else if(vseed == 21){
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
    else if(vseed == 87){
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
    else if(vseed == 88){
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
    else if(vseed == 89){
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
	// ***** For BPI vseed = 114 and 116 with [110]; and for BPII vseed = 115 and 119 with [111] ***** //////////
	// ***** vseeds 124, 126, 125, 129 and so on are combinations (with a 1 at the beginning) of vseed 2 o 3 (cholesteric)
	// ***** and BPI's vseeds 4 & 6, and BPII's vseeds 5 & 9.
	// *** for vseed 141, 142, 143 BPII 100, 110, 111 outter shells with BPIII random inner core //

	else if(vseed == 114 || vseed == 116 || vseed == 115 || vseed == 119 ||	
			vseed == 124 || vseed == 126 || vseed == 125 || vseed == 129 ||
			vseed == 134 || vseed == 136 || vseed == 135 || vseed == 139 || 
			vseed == 141 || vseed == 142 || vseed == 143 || 
			vseed == 874 || vseed == 875 || vseed == 876 || vseed == 879)
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

		if( vseed == 114 || vseed == 116 || vseed == 124 || vseed == 126 || vseed == 134 || vseed == 136 || 
			vseed == 141 || vseed == 142 || vseed == 143 ||
			vseed == 874 || vseed == 876){
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
							if(vseed == 114 || vseed == 116 || vseed == 115 || vseed == 119){
								
								
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
							else if(vseed == 141 || vseed == 142 || vseed == 143){
								
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
							else if(vseed == 124 || vseed == 126 || vseed == 125 || vseed == 129 ||
									vseed == 134 || vseed == 136 || vseed == 135 || vseed == 139){
					
								x = (i-rx)*dx;
								y = (j-ry)*dy;
								z = (k-rz)*dz;
								dis = sqrt(x*x+y*y+z*z);

								if(vseed == 124 || vseed == 126 || vseed == 125 || vseed == 129){
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
							else if(vseed == 874 || vseed == 875 || vseed == 876 || vseed == 879){
								
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
							if(vseed == 114 || vseed == 124 || vseed == 134 || vseed == 874){

								x = (i - Nx * 0.5) * cst * isq2;
								y = (j - Ny * 0.5) * cst * isq2;
								z = (k - Nz * 0.5) * cst * isq2;
						
								Qold[nd * 6 + 0] = A * (- sin(y) * cos(x) - sin(x) * cos(z) + 2 * sin(z) * cos(y));
								Qold[nd * 6 + 3] = A * (- sin(z) * cos(y) - sin(y) * cos(x) + 2 * sin(x) * cos(z));
								Qold[nd * 6 + 5] = A * (- sin(x) * cos(z) - sin(z) * cos(y) + 2 * sin(y) * cos(x));
								Qold[nd * 6 + 1] = A * (- sq2 * sin(x) * sin(z) - sq2 * cos(y) * cos(z) + sin(x) * cos(y));
								Qold[nd * 6 + 2] = A * (- sq2 * sin(z) * sin(y) - sq2 * cos(x) * cos(y) + sin(z) * cos(x));
								Qold[nd * 6 + 4] = A * (- sq2 * sin(y) * sin(x) - sq2 * cos(z) * cos(x) + sin(y) * cos(z));
							}

							else if(vseed == 115 || vseed == 125 || vseed == 135 || vseed == 141 || vseed == 875){
										
								x = i - Nx * 0.5;
								y = j - Ny * 0.5;
								z = k - Nz * 0.5;

								Qold[nd * 6 + 0] = A * (cos(cst * z) - cos(cst * y));
								Qold[nd * 6 + 1] = A * sin(cst * z);
								Qold[nd * 6 + 2] = A * sin(cst * y);
								Qold[nd * 6 + 3] = A * (cos(cst * x) - cos(cst * z));
								Qold[nd * 6 + 4] = A * sin(cst * x);
								Qold[nd * 6 + 5] = A * (cos(cst * y) - cos(cst * x));
							}

							else if(vseed == 116  || vseed == 126 || vseed == 136 || vseed == 876){

								xi = (i - Nx * 0.5) * cst * isq2;
								yi = (j - Ny * 0.5) * cst * isq2;
								zi = (k - Nz * 0.5) * cst * isq2;
								
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
													
							else if (vseed == 142){
							
								xi = i - Nx * 0.5;
                        		yi = j - Ny * 0.5;
                            	zi = k - Nz * 0.5;

								x = xi;
								y = cos(theta) * yi + sin(theta) * zi;
								z = -sin(theta) * yi + cos(theta) * zi;

								Qold[nd * 6 + 0] = A * (cos(cst * z) - cos(cst * y));
								Qold[nd * 6 + 1] = A * sin(cst * z);
								Qold[nd * 6 + 2] = A * sin(cst * y);
								Qold[nd * 6 + 3] = A * (cos(cst * x) - cos(cst * z));
								Qold[nd * 6 + 4] = A * sin(cst * x);
								Qold[nd * 6 + 5] = A * (cos(cst * y) - cos(cst * x));}

							else if(vseed == 119  || vseed == 129 || vseed == 139 || vseed == 143 || vseed == 879){

								xi = i - Nx * 0.5;
								yi = j - Ny * 0.5;
								zi = k - Nz * 0.5;
													
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
