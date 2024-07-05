#include"finite.h"

int main(int argc, char * argv[]){

	//read_param();
		
        MPI_Init(&argc, &argv);
        MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &shmcomm);
        MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
        MPI_Comm_rank(MPI_COMM_WORLD, &myid);

	int i;

	flag = true;
	flag_2 = false;

	double time_spend;
	double begin, end;
	double deltat;
	double shell_accuracy = 10e-1;

	FILE* energy;

	begin = MPI_Wtime();
	//read in parameters
	

	if(!read_param()){
			MPI_Comm_free(&shmcomm);
			int MPI_Finalize();
			return 0;
	}

	dt = tmin;
	dE = 1;
	el_old = 1;
	cycle = 0;
	deltat = (tmax - tmin) / increment;

	S = 0.25 * (1.0 + 3.0 * sqrt(1.0 - 8.0 / (3.0 * U)));
	//	printf("Theoretical value is %lf.\n", third * (1 - third * U) * S * S - 2 * third * third * third * S * S * S * U + U / 9 * S * S * S * S );
	
	S2 = 0.25 * (1.0 + 3.0 * sqrt(1.0 - 8.0 / (3.0 * U2)));

	if(myid == root){
		//define droplet and boundary, introduce particles, initialize qtensor;
		if(!initial()){
			flag = false;
			return 0;
		}		
	}
	
	MPI_Bcast(&flag, 1, MPI_BYTE, root, MPI_COMM_WORLD);	
	
	
	//For all infomation initialized in root processor, scatter them to all other processors.
	if(!scatter()){
		flag = false;
	}
	
	if(seed == -1442 || seed == -1443){
		flag_2 = true;
	}
		
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Win_fence(0, win);
	
	//Evolution
	while(flag){
		//Every 1000 steps calculate energies.
		if(cycle % check_every == 0){
			free_energy();

			if(geo == -44){

				if(fabs(dE) < shell_accuracy && flag_2 == false && (seed != -1442 || seed != -1443)){
					flag_2 = true;
				}
				else if(fabs(dE) < accuracy){
					flag = false;
				}
			}
			else{
				if(fabs(dE) < accuracy){
					flag = false;
				}
			}
			
		}

		if(cycle % trace_checker == 0){ 
			//Every 10000 steps check the trace of Qtensor
			if(myid == root){	
				for(i = 0; i < droplet; i++){
					//				checktr(&q[i * 6]);
					if(!checktr(&q[i * 6])){
						flag = false;
						printf("%d\n", i);
					}
				}
				if(!flag){
					printf("Error in the trace of q; cycle : %d.\n", cycle);
					//MPI_Bcast(&flag, 1, MPI_BYTE, root, MPI_COMM_WORLD);
				}
				//output();
			}	
			MPI_Bcast(&flag, 1, MPI_BYTE, root, MPI_COMM_WORLD);	
		}

		//Mi modificación para evitar que escriba cada mil ciclos en el disco duro.
		if(cycle % save_every == 0 && myid == root){
			if(geo != -44){
				output();
			}
			else{
				evolving_output(cycle);
			}
		}

		//Wait until all the processors are ready and relax the system, first bulk and then boundary.
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Win_fence(0, win);
		if(flag){
			if(geo == -44 && flag_2 == true){
				relax_bulk();
			}
			else if(geo == -44 && flag_2 == false){
				relax_evolving_bulk();	
			}
			else{
				relax_bulk();
			}
		} 

		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Win_fence(0, win);
		if(flag){
			if(geo == -44 && flag_2 == true){
				//relax_evolving_surf();
			}
			else{
				relax_surf();
			}
		} 
		
/*
		if(myid == root){	
			printf("check2.\n");
		}
*/
		if(dt < tmax){
			dt += deltat;
			if(dt >= tmax)	dt = tmax;
		}

		cycle++;
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Win_fence(0, win);
		
	}

	free_energy();

	end = MPI_Wtime();

	if(myid == root){
		output();

		//Calculate time used.
		time_spend = (double)(end - begin) / 60.0;

		energy = fopen("energy.out", "a");
		if(time_spend < 60){
			fprintf(energy, "\nTime used:	%lf min.\n", time_spend);
			printf("\nTime used:	%lf min.\n", time_spend);
		}
		else{
			fprintf(energy, "\nTime used:	%lf h.\n", time_spend / 60.0);
			printf("\nTime used:	%lf h.\n", time_spend / 60.0);
		}
		fclose(energy);	
	}

	//deallocate dynamic arrays
	free_q();
	MPI_Win_free(&win);
	MPI_Win_free(&win2);
    MPI_Comm_free(&shmcomm);
    MPI_Finalize();

	return 0;
}
