#include "read_param.h"

#define READ_PARAM_LINE(expected_count, ...)            \
	do {                                                \
		char line[512];                                 \
		if (fgets(line, sizeof(line), param) == NULL) { \
			printf("Error reading param.in.\n");        \
			fclose(param);                              \
			return false;                               \
		}                                               \
		if (sscanf(line, __VA_ARGS__) != expected_count) { \
			printf("Error reading param.in.\n");        \
			fclose(param);                              \
			return false;                               \
		}                                               \
	} while (0)

bool read_param() {

	//seed = 0;
	FILE* param = fopen("param.in", "r");
	int degenerate_in = 0;
	int uppersurf_in = 0;
	int lowersurf_in = 0;
	int surfdegen_in = 0;
	int DoubleU_in = 0;

	if (param == (FILE*)NULL) {
		printf("No param.in file found!\n");
		return false;
	}

	READ_PARAM_LINE(1, "Nx %d", &Nx);
	READ_PARAM_LINE(1, "Ny %d", &Ny);
	READ_PARAM_LINE(1, "Nz %d", &Nz);
	READ_PARAM_LINE(1, "Lx %lf", &Lx);
	READ_PARAM_LINE(1, "Ly %lf", &Ly);
	READ_PARAM_LINE(1, "Lz %lf", &Lz);
	READ_PARAM_LINE(1, "W %lf", &W);
	READ_PARAM_LINE(1, "U %lf", &U);
	READ_PARAM_LINE(1, "L1 %lf", &L1);
	READ_PARAM_LINE(1, "L2 %lf", &L2);
	READ_PARAM_LINE(1, "L3 %lf", &L3);
	READ_PARAM_LINE(1, "L4 %lf", &L4);
	READ_PARAM_LINE(1, "chiral %d", &chiral);
	READ_PARAM_LINE(1, "qch %lf", &qch);
	READ_PARAM_LINE(1, "redshift %lf", &redshift);
	READ_PARAM_LINE(1, "geo %d", &geo);
	READ_PARAM_LINE(1, "degenerate %d", &degenerate_in);
	READ_PARAM_LINE(1, "tiltAngle %lf", &tiltAngle);
	READ_PARAM_LINE(1, "infinite %d", &infinite);
	READ_PARAM_LINE(1, "Np %d", &Np);
	READ_PARAM_LINE(1, "Rp %lf", &Rp);
	READ_PARAM_LINE(1, "Wp %lf", &Wp);
	READ_PARAM_LINE(1, "seed %d", &seed);
	READ_PARAM_LINE(1, "rand_seed %d", &rand_seed);
	READ_PARAM_LINE(2, "tmin tmax %lf %lf", &tmin, &tmax);
	READ_PARAM_LINE(1, "increment %lf", &increment);
	READ_PARAM_LINE(1, "accuracy %lf", &accuracy);
	READ_PARAM_LINE(3, "init_dir %lf %lf %lf", &init_dir[0], &init_dir[1], &init_dir[2]);
	READ_PARAM_LINE(3, "dir1 %lf %lf %lf", &dir1[0], &dir1[1], &dir1[2]);
	READ_PARAM_LINE(3, "dir2 %lf %lf %lf", &dir2[0], &dir2[1], &dir2[2]);
	READ_PARAM_LINE(1, "UpperSurface %d", &uppersurf_in);
	READ_PARAM_LINE(1, "LowerSurface %d", &lowersurf_in);
	READ_PARAM_LINE(1, "LowerSurfaceDegen %d", &surfdegen_in);
	READ_PARAM_LINE(1, "DoubleU Mode %d", &DoubleU_in);
	READ_PARAM_LINE(1, "U2 %lf", &U2);
	READ_PARAM_LINE(1, "iRx %lf", &iRx);
	READ_PARAM_LINE(1, "iRy %lf", &iRy);
	READ_PARAM_LINE(1, "iRz %lf", &iRz);
	READ_PARAM_LINE(1, "Save Every %d", &save_every);
	READ_PARAM_LINE(1, "Check Every %d", &check_every);
	READ_PARAM_LINE(1, "Stop At %d", &stopat);
	READ_PARAM_LINE(1, "Check Trace At %d", &trace_checker);

		degenerate = (degenerate_in != 0);
		uppersurf = (uppersurf_in != 0);
		lowersurf = (lowersurf_in != 0);
		surfdegen = (surfdegen_in != 0);
		DoubleU = (DoubleU_in != 0);
		if(myid == root){
		printf("Nx %d\n", Nx);
		printf("Ny %d\n", Ny);
		printf("Nz %d\n", Nz);
		printf("Lx %.1lf\n", Lx);
		printf("Ly %.1lf\n", Ly);
		printf("Lz %.1lf\n", Lz);
		printf("W %lf\n", W);
		printf("U %lf\n", U);
		printf("L1 %.3lf\n", L1);
		printf("L2 %.3lf\n", L2);
		printf("L3 %.3lf\n", L3);
		printf("L4 %.3lf\n", L4);
		printf("chiral %d\n", chiral);
		printf("qch %lf\n", qch);
		printf("redshift %lf\n", redshift);
		printf("geo %d\n", geo);
		printf("degenerate %d\n", degenerate);
		printf("tilt angle = %.2lf\n", tiltAngle);
		printf("infinite %d\n", infinite);
		printf("Np %d\n", Np);
		printf("Rp %lf\n", Rp);
		printf("Wp %lf\n", Wp);
		//	printf("pdegenerate %d\n", pdegenerate);
		//	printf("pinfinite %d\n", pinfinite);
		printf("seed %d\n", seed);
		printf("rand_seed %d\n", rand_seed);
		printf("tmin tmax %.e %.e\n", tmin, tmax);
		printf("increment %.e\n", increment);
		printf("accuracy %.e\n", accuracy);
		printf("init_dir is %.1lf, %.1lf, %.1lf\n", init_dir[0], init_dir[1], init_dir[2]);
		printf("dir1 is %.1lf, %.1lf, %.1lf\n", dir1[0], dir1[1], dir1[2]);
		printf("dir2 is %.1lf, %.1lf, %.1lf\n", dir2[0], dir2[1], dir2[2]);
		printf("Upper Surface is %d\n", uppersurf);
		printf("Lower Surface is %d\n", lowersurf);
		printf("Checkpoint every %d!\n", save_every);
		printf("Energy will be compared every %d!\n", check_every);
		printf("Trace will be checked every %d!\n", trace_checker);
		if(stopat != 0){
			printf("Job will be\033[1;31m STOPPED\033[0m after %d!\n", stopat);
		}
		else{
			printf("Job running until dE is reached!\n");
		}

		if (surfdegen) {
			printf("Lower Surface Degenerate will degenerate\n");
		}
		else {
			printf("Lower Surface Degenerate will NOT degenerate\n");
		}		

		if (infinite && degenerate) {
			printf("Degenerate planar anchoring cannot be infinite.\n");
			return false;
		}

		if (!norm_v(dir1) || !norm_v(dir2)) {
			printf("Invalid dir1 or dir2 input.\n");
			return false;
		}

		if (tmax < tmin) {
			printf("Error for time input.\n");
			return false;
		}

		if (!norm_v(init_dir)) {
			printf("Problems in initial direction!\n");
			return false;
		}

		if (lowersurf == 0 && dir2[2] != 0) {
			printf("Error in dir2! Planar can only be in x or z direction.\n");
			return false;
		}
		/*Por ahora dejaré el comando así, pero después podremos configurarlo para que pueda
		configuraciones inclinadas.		*/
		else if (lowersurf == 1 && (dir2[0] != 0 || dir2[1] != 0)) {
			printf("Error in dir2! Homeotropic can only be in z direction.\n");
			return false;
		}

		//por ahora mantendremos normales los vectores en la superficie superior.
		if (uppersurf == 1 && (dir1[0] != 0 && dir1[1] != 0)) {
			printf("Error in dir1! Homeotropic can't be in x or z direction.\n");
			return false;
		}
		else if (uppersurf == 0 && dir1[2] == 1) {
			printf("Error in dir1! Planar can't be in y direction.\n");
			return false;
		}

		if (DoubleU) {
			printf("Double U mode is activated.\n");
			printf("U for outer shell = %lf\n", U2);
			printf("Inner Radii iRx = %lf; iRy = %lf; iRz = %lf\n", iRx, iRy, iRz);
		}
		else {
			printf("Double U mode is NOT activated.\n");
		}
	}
	fclose(param);
	return true;
	
}

bool read_nppos(double **pos){
	int Npnum, i, j;
	int coor;
	double x1, y1, z1, x2, y2, z2;
        FILE *nppos;
	nppos = fopen("nppos.in","r");
	if(nppos == (FILE*)NULL){
		printf("File nppos.in not found.\n");
		return false;
	}
	//read in NP number Np 
	//read in Coordinate: 0 for cartesian; 1 for spherical
        fscanf(nppos, "Np\t%d\nCoor\t%d\n", &Npnum, &coor);
	if(Npnum != Np){
		printf("Error in nanoparticle position file (NP number %d and %d conflict.)\n", Np, Npnum);
		return false;
	}
	if(coor != 0 && coor != 1){
		printf("Invalid input for nppos.in.\nSpherical coordinate: input 1; Cartesian coordinate: input 0.\n");
		return false;
	}
	fscanf(nppos, "posx\tposy\tposz\tanchoring\n");
	//read in particle position and anchoring(pos[i][3]), if spherical coordinate, change to cartesian
	//anchoring: 
	//0 for infinite anchoring
	//1 for noninfinite homeotropic
	//2 for planar degenerate
	for(i = 0; i < Np; i ++){
		fscanf(nppos, "%lf\t%lf\t%lf\t%lf\n", &pos[i][0], &pos[i][1], &pos[i][2], &pos[i][3]);
		if(pos[i][3] > 2 || pos[i][3] < 0){
			printf("Wrong input for anchoring condition.\n");
			return false;
		}
		if(coor == 1){
			x1 = lrint(Nx * 0.5) + pos[i][0] * sin(pos[i][1] / 180.0 * M_PI) * cos(pos[i][2] / 180.0 * M_PI);
			y1 = lrint(Ny * 0.5) + pos[i][0] * sin(pos[i][1] / 180.0 * M_PI) * sin(pos[i][2] / 180.0 * M_PI);
			z1 = lrint(Nz * 0.5) + pos[i][0] * cos(pos[i][1] / 180.0 * M_PI);
			pos[i][0] = x1;
			pos[i][1] = y1;
			pos[i][2] = z1;
		}
/*		if((pos[i][2] + Rp > Nz - 2) || (pos[i][2] - Rp < 2)){
			printf("Particle colloids with the wall.\n");
			fclose(nppos);
			return false;
		}
		if((pos[i][1] + Rp > Ny - 2) || (pos[i][1] - Rp < 2)){
			printf("Particle colloids with the wall.\n");
			fclose(nppos);
			return false;
		}
		if((pos[i][0] + Rp > Nx - 2) || (pos[i][0] - Rp < 2)){
			printf("Particle colloids with the wall.\n");
			fclose(nppos);
			return false;
		}
*/
		printf("%lf\t%lf\t%lf\t%lf\n", pos[i][0], pos[i][1], pos[i][2], pos[i][3]);
	}
	//examine if particles overlap
	if(Np > 1){
		for(i = 0; i < Np - 1; i ++){
			for(j = i + 1; j < Np; j ++){
				x1 = pos[i][0] - pos[j][0];
				y1 = pos[i][1] - pos[j][1];
				z1 = pos[i][2] - pos[j][2];
				if((x1 * x1) + (y1 * y1) + (z1 * z1) <= (4 * (Rp + 2) * (Rp + 2))){
					printf("Particle %d and %d overlap.\n", i, j);
					printf("Distance is %lf.\n\n", sqrt(x1 * x1 + y1 * y1 + z1 * z1));
					return false;
				}
			}
		}
	}	

	//print nanoparticle positions
	fclose(nppos);
	return true;
}
