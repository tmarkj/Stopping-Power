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

#endif
