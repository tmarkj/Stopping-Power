#ifndef STOPPING_POWER_H
#define STOPPING_POWER_H

#define TABLE_LENGTH 1000
static const unsigned int HEADER_SIZE = 4;

struct Stopping_power {

	char* ion;
	char* filter_material;
	double energy_array[TABLE_LENGTH];
	double stop_pow[TABLE_LENGTH];
	double range[TABLE_LENGTH];

};


int read_table(char* ion, char* filter_material, double* energy, double* stop_pow);

double interp1d(double * x_array, double* y_array, double x_point, int length);

int initialize_stopping_power(struct Stopping_power * stop_pow_struct);

void range_cumtrapz(double* range, double* energy_array, double* stop_pow_array, int length);

double E_out(struct Stopping_power * stop_pow_struct, double E_in, double thickness);

double E_in(struct Stopping_power * stop_pow_struct, double E_out, double thickness);

void get_centered_from_edge(double * edge_array, double * centered_array, int edge_len);

void get_edge_from_centered(double * centered_array, double * edge_array, int center_len);

void diff_stop(double * in_array, double * out_array, int in_len);

void E_out_spectrum(struct Stopping_power * stop_pow_struct, double * E_in_array, double * E_out_array,
	double * yield_in_array, double * yield_out_array, double thickness, int len);

#endif
