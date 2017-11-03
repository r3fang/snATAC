// gcc -g -O2 calculate_jacard_index_matrix.c -o calculate_jacard_index_matrix
// wc -l in.mat # num_row # num_row = 12932
// awk '{print NF}' in.mat | sort -nu | tail -n 1 # num_col = 140102
#include <stdio.h>
#include <string.h>

#ifndef PACKAGE_VERSION
#define PACKAGE_VERSION "0.7.23-r15"
#endif

int main(int argc, char *argv[])
{
	int i, j, k, num_row, num_col;	
	if (argc == 1) {  
	        fprintf(stderr, "Usage: %s <in.mat> <num_row> <num_col> <out.txt> \n", argv[0]);  
	        return 1;  
	}  
	num_row = strtol(argv[2], NULL, 10);
	num_col = strtol(argv[3], NULL, 10);

	int** mat=malloc(num_row*sizeof(int*)); 
	for(i=0;i<num_row;++i)
	mat[i]=malloc(num_col*sizeof(int));
	int num; 
	FILE *file;
	file=fopen(argv[1], "r");
	for(i = 0; i < num_row; i++){
		for(j = 0; j < num_col; j++){
			if (!fscanf(file, "%d", &mat[i][j])) 
				break;
	      }
	  }
	 fclose(file);

 	double** jacard=malloc(num_row*sizeof(double*)); 
	double union_ij;
	double joint_ij;
	for(i=0;i<num_row;++i)
	jacard[i]=malloc(num_row*sizeof(double));
 	for(i = 0; i < num_row; i++){
		for(j = i; j < num_row; j++){
			if(i == j){
				jacard[i][j] = 0;
			}else{
				union_ij = 0;
				joint_ij = 0;
				for(k = 0; k < num_col; k ++){
					if(mat[i][k] + mat[j][k] > 0) union_ij += 1;
					if(mat[i][k] + mat[j][k] == 2) joint_ij += 1;
				}
				jacard[i][j] = (joint_ij+1)/(union_ij+1);
				jacard[j][i] = (joint_ij+1)/(union_ij+1);
			}
		}
	}	
	
	FILE *fout;
	fout=fopen(argv[4], "w");
 	for(i = 0; i < num_row; i++){
		for(j = 0; j < num_row; j++){
 			fprintf(fout, "%.2f\t", jacard[i][j]*1000);			
		}
		fprintf(fout, "\n");
	}
	
	fclose(fout);  
}