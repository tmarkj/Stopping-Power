import numpy as np
from scipy.integrate import cumtrapz
from scipy import interpolate
import os

class Stopping_power:
    """
    Class for doing stopping power calculations.
    """

    def __init__(self, ion, filter_material):
        """
        Parameters 
        ----------
        ion : str
            The type of particle that will pass through the material.
            Ex: "H" for hydrogen, "D" for deuteron, "T" for triton

        filter_material : str
            The material of that the high energy ions will pass through
            Ex: "Ta" for tantalum, "Al" for aluminum
        """

        self.ion = ion
        self.filter_material = filter_material

        # The SRIM output tables are stored in the Tables/ directory
        module_dir = os.path.dirname(os.path.realpath(__file__))
        SRIM_file = f"Tables/{ion}_in_{filter_material}"
        raw_SRIM_data = np.genfromtxt(f"{module_dir}/{SRIM_file}", skip_header = 4) 
        
        self.ion_E = np.asarray(raw_SRIM_data[:,0], dtype=float)/1e3 # MeV
        elec_dEdx = np.asarray(raw_SRIM_data[:,1], dtype=float)
        nuc_dEdx = np.asarray(raw_SRIM_data[:,2], dtype=float)
        net_dEdx = elec_dEdx + nuc_dEdx

        # Integrate the stopping power to find the range of the particle
        # as a function of energy
        self.range_array = cumtrapz(1/net_dEdx, self.ion_E, initial=0.0)*1e3 #um

        self.interp_range = interpolate.interp1d(self.ion_E, self.range_array)
        self.interp_energy = interpolate.interp1d(self.range_array, self.ion_E)

    def E_out(self, E_in, thickness):
        """
        Finds the energy of the ion after it passes through the material.
        It will give an error message and return np.nan of the filter thickness
        exceeds the range of the particle in the material.

        Parameters
        ----------
        E_in : float
            The energy of the ion before the filter. This is in units of MeV.

        thickness : float
            This is the thickness of the filter material in um.

        Returns
        -------
        E_out: float (or nan if it's ranged out)
            The energy of the particle after it passes through the material.
            Units of MeV.
        """
        if thickness >= self.range(E_in):
            print(f"Range of {E_in:.3f} MeV {self.ion} in {self.filter_material} is less than the thickness! Particle is ranged out.")
            return np.nan
        else:
            return self.interp_energy(self.interp_range(E_in) - thickness)

    def E_in(self, E_out, thickness):
        """
        Finds the energy that the ion had before it passed through the
        material.

        Parameters
        ----------
        E_out : float
            The energy of the ion after it has passed through the material.
            Units are MeV.

        thickness : float
            This is the thickness of the filter material in um.

        Returns
        -------
        E_in : float
            The energy of the particle before it passed through the material.
            Units are MeV.
        """
        return self.interp_energy(thickness + self.interp_range(E_out))

    def range(self, E):
        """
        Finds the range of the ion in the material at a given energy.

        Parameters
        ----------
        E : float
            The energy of the ion in MeV

        Returns
        -------
        range : float
            The range of the ion at the given energy in the material. Units are
            um.
        """
        return self.interp_range(E)
