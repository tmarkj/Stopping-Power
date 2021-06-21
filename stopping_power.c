#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "stopping_power.h"



int read_table(char* ion, char* filter_material, double* energy, double* stop_pow) {

	char buf[128];
	char filepath[128];
	int counter = 0;
	char *energy_in, *nuc_pow_in, *ele_pow_in;
	const char* delimiter = "   ";

	snprintf(filepath, sizeof(filepath), "%s/Tables/%s_in_%s", TABLE_PATH, ion, filter_material);

	FILE * file = fopen(filepath, "r");
	
	if (file == NULL){
		fprintf(stderr, "file %s cannot open!\n", filepath);
		return -1;
	}
	
	while(!feof(file)){
		if(fgets(buf, 127, file) != NULL) {
			if (counter < HEADER_SIZE){
				++counter;
			} else {
				energy_in = strtok(buf, delimiter);
				ele_pow_in = strtok(NULL, delimiter);		
				nuc_pow_in = strtok(NULL, delimiter);		
				
				*energy++ = atof(energy_in)/1000.0; // convert from keV to MeV
				*stop_pow++ = atof(ele_pow_in) + atof(nuc_pow_in);

				++counter;
			}
		}
	}
	
	fclose(file);
		
	return 1;
}

void range_cumtrapz(double* range, double* energy_array, double* stop_pow_array, int length){
	double diff;

	range[0] = 0;
	
	for (int i = 1; i < length; ++i) {
		diff = energy_array[i] - energy_array[i-1];
		range[i] = 1/2.0*(1/stop_pow_array[i] + 1/stop_pow_array[i-1])*diff*1e3 + range[i-1];

	}
}




double interp1d(double * x_array, double* y_array, double x_point, int length){
	int close_index; // and also smaller than x

	// Check to see if value is within the array
	// This assumes that x_array is monotonically increasing (as it should be!)
	// Could use a binary search to make this faster but I'm incompetent
	if (x_point < x_array[0] || x_point > x_array[length-1]){
		printf("Interpolation error! Outside the range of the x array.\n");
		printf("min %f \t max %f \t point %f \n", x_array[0], x_array[length], x_point);
		return -1;
	}

	for (int i = 0; i < length; ++i){
		if (x_point > x_array[i]){
			close_index = i;
		} else {
			break;
		}
	}

	double slope = (y_array[close_index+1] - y_array[close_index])/(x_array[close_index+1] - x_array[close_index]);

	return y_array[close_index] + (x_point - x_array[close_index])*slope;
	

}



int initialize_stopping_power(struct Stopping_power * stop_pow_struct){

	read_table(stop_pow_struct->ion, stop_pow_struct->filter_material, stop_pow_struct->energy_array, stop_pow_struct->stop_pow);

	return 1;	

}

double E_out(struct Stopping_power * stop_pow_struct, double E_in, double thickness){

	double range_at_E_in = interp1d(stop_pow_struct->energy_array, stop_pow_struct->range, E_in, TABLE_LENGTH);

	if (range_at_E_in < thickness){
		return 0;
	}

	return interp1d(stop_pow_struct->range, stop_pow_struct->energy_array, range_at_E_in - thickness, TABLE_LENGTH);
}

double E_in(struct Stopping_power * stop_pow_struct, double E_out, double thickness){
	double range_at_E_out = interp1d(stop_pow_struct->energy_array, stop_pow_struct->range, E_out, TABLE_LENGTH);

	return interp1d(stop_pow_struct->range, stop_pow_struct->energy_array, range_at_E_out + thickness, TABLE_LENGTH);
}


void get_centered_from_edge(double * edge_array, double * centered_array, int edge_len){

	double dx = edge_array[1] - edge_array[0];

	for (int i = 0; i < edge_len - 1; ++i) {
		centered_array[i] = edge_array[i] + dx/2.0;
	}
}

void get_edge_from_centered(double * centered_array, double * edge_array, int center_len){
	
	double dx = centered_array[1] - centered_array[0];
	edge_array[0] = centered_array[0] - dx/2.0;

	for (int i = 0; i < center_len; ++i){
		edge_array[i+1] = centered_array[i] + dx/2.0;
	}

}

void diff_stop(double * in_array, double * out_array, int in_len){
	// Note that lengths of the arrays are different!
	// in_array is len n, but out_array is len n-1 
	// It also has the property that for the first non-zero diff,
	// it will make that a large number.
	// This is purely for stopping power weirdness purposes
	// so don't use this for anything else

	for (int i = 0; i < in_len - 1; ++i){
		if ((in_array[i+1] != 0.0) && (in_array[i] == 0.0) ) {
			out_array[i] = 1e10;
		} else {
			out_array[i] = in_array[i+1] - in_array[i];		
		}
	}

}


void E_out_spectrum(struct Stopping_power * stop_pow_struct, double * E_in_array, double * E_out_array, 
							double * yield_in_array, double * yield_out_array, double thickness, int len){	

	double E_out_value;
	int index_first_valid_energy;

	double * E_in_array_edge = (double *) malloc((len+1)*sizeof(double));
	get_edge_from_centered(E_in_array, E_in_array_edge, len);

	double * E_out_array_edge = (double *) malloc((len+1)*sizeof(double));
	
	for (int i = 0; i < len + 1; ++i) {
		E_out_value = E_out(stop_pow_struct, E_in_array_edge[i], thickness);

		// Check to see if the E_out is valid
		if (E_out_value < .005) E_out_value = 0;
			
		E_out_array_edge[i] = E_out_value; 
	}

	double * dE_in_array = (double *)malloc(len*sizeof(double));
	double * dE_out_array = (double *)malloc(len*sizeof(double));

	diff_stop(E_in_array_edge, dE_in_array, len+1);
	diff_stop(E_out_array_edge, dE_out_array, len+1);

	for (int i = 0; i < len; ++i){
		if (dE_out_array[i] == 0.0) {
			yield_out_array[i] = 0.0;
		} else{
			yield_out_array[i] = yield_in_array[i]*dE_in_array[i]/dE_out_array[i];
		}
	}

	get_centered_from_edge(E_out_array_edge, E_out_array, len+1);


	free(E_in_array_edge);
	free(E_out_array_edge);
	free(dE_in_array);
	free(dE_out_array);


}


