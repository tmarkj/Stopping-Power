import numpy as np
from scipy.integrate import cumtrapz
from scipy import interpolate

class Stopping_power:

    def __init__(self, ion, filter_material):

        self.ion = ion
        self.filter_material = filter_material

        SRIM_file = f"Tables/{ion}_in_{filter_material}"
        raw_SRIM_data = np.genfromtxt(SRIM_file, skip_header = 4) 
        
        self.ion_E = np.asarray(raw_SRIM_data[:,0], dtype=float)/1e3 # MeV
        elec_dEdx = np.asarray(raw_SRIM_data[:,1], dtype=float)
        nuc_dEdx = np.asarray(raw_SRIM_data[:,2], dtype=float)
        net_dEdx = elec_dEdx + nuc_dEdx

        self.range_array = cumtrapz(1/net_dEdx, self.ion_E, initial=0.0)*1e3 #um

        self.interp_range = interpolate.interp1d(self.ion_E, self.range_array)
        self.interp_energy = interpolate.interp1d(self.range_array, self.ion_E)

    def E_out(self, E_in, thickness):
        if thickness >= self.range(E_in):
            print(f"Range of {E_in:.3f} MeV {self.ion} in {self.filter_material} is less than the thickness! Particle is ranged out.")
            return np.nan
        else:
            return self.interp_energy(self.interp_range(E_in) - thickness)

    def E_in(self, E_out, thickness):
        return self.interp_energy(thickness + self.interp_range(E_out))

    def range(self, E):
        return self.interp_range(E)
