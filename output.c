#include "finite.h"

void output(){
	FILE* energy;
	FILE* Qtensor;
	int i, j, k, l, indx, n;
	double zero = 0;

	//print energy
	energy = fopen("energy.out", "a");
	if(DoubleU){
		fprintf(energy,"%d\t%.9lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", cycle, dE, en_ldg[0], en_ldg[1], en_ldg[2], en_el[0], en_el_in[0], en_el_out[0], en_el[4], en_el_in[4], en_el_out[4], en_surf[0], en_tot);
	}
	else{
		fprintf(energy,"%d\t%.9lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", cycle, dE, en_ldg[0], en_el[0], en_el[1], en_el[2], en_el[3], en_el[4], en_surf[0], en_surf[1], en_tot);
	}
	fclose(energy);

	//print Qtensor
	Qtensor = fopen("Qtensor.bin", "wb");
	indx = 0;
	for(l = 0; l < tot; l++){
		if(drop[l] || boundary[l] || nboundary[l]){
			fwrite(&q[indx * 6], sizeof(double), 6, Qtensor);
			indx ++;
		}
	}
	fclose(Qtensor);
}
