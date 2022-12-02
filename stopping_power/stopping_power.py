import numpy as np
from scipy.integrate import cumtrapz
from scipy import interpolate
import os

    
def get_centered_from_edge(edge_array):
    """
    dx = edge_array[1] - edge_array[0]
    centered_array = []

    for i in range(len(edge_array) - 1):
        centered_array.append(edge_array[i] + dx/2)
    """
    dx = np.gradient(edge_array)
    centered_array = []

    for i in range(len(edge_array) - 1):
        centered_array.append(edge_array[i] + dx[i]/2)

    return np.asarray(centered_array)


def get_edges_from_centered(center_array):
    """
    dx = center_array[1] - center_array[0]
    edge_array = [center_array[0] - dx/2]

    for i in range(len(center_array)):
        edge_array.append(center_array[i] + dx/2)
    """

    dx = np.gradient(center_array)
    edge_array = [center_array[0] - dx[0]/2]

    for i in range(len(center_array)):
        edge_array.append(center_array[i] + dx[i]/2)

    return np.array(edge_array)


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

        self.interp_range = interpolate.interp1d(self.ion_E, self.range_array, bounds_error=False)
        self.interp_energy = interpolate.interp1d(self.range_array, self.ion_E, bounds_error=False)

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

        ## So the code can handle array an non-array input
        not_array = False
        E_in_array = np.array(E_in)
        if E_in_array.shape == ():
            not_array = True
            E_in_array = np.array([E_in])

        ## Keep track of which parts get ranged out and give NaN
        ranged_out_indices = np.argwhere(thickness >= self.range(E_in_array))
        good_indices = np.where(thickness < self.range(E_in_array))
        if len(ranged_out_indices) > 0:
            #print(f"Range of {self.ion} in {self.filter_material} is less than the thickness! Particle is ranged out.")
            pass


        E_out = np.zeros(E_in_array.shape)
        E_out[good_indices] = self.interp_energy(self.interp_range(E_in_array[good_indices]) - thickness)
        E_out[ranged_out_indices] = np.nan

        ## Return the same type as the input (array or not array)
        if not_array:
            return E_out[0]
        else:
            return E_out

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

    def thickness(self, E_in, E_out):
        """
        Finds the thickness that took particle from energy E_in to energy E_out 

        Parameters
        ----------
        E_in : float
            The energy of the ion into the material in MeV

        E_out : float
            The energy of the ion out of the material in MeV

        Returns
        -------
        thickness : float
            The thickness of material required for take energy E_in to energy E_out.
            Units are um.
        """
        return self.interp_range(E_in) - self.interp_range(E_out)


    def E_out_spectrum(self, E_in_array, yields_in_array, thickness, errorbars=None):
        """
        Finds the spectrum after the particles passed through the material.

        Parameters
        ----------
        E_in_array : array of floats (numpy array)
            The energy array for the spectrum. MeV units.

        yields_in_array : array of floats (numpy array)
            The array of yields per MeV of the in spectrum.

        Returns
        -------
        E_out_array : array of floats (numpy array)
            The energy of the out spectrum. MeV units.

        yields_out_array : array of floats (numpy array)
            The array of yields per MeV of the out spectrum .

        """

        E_in_array_edges = get_edges_from_centered(E_in_array)
        E_out_array_edges = self.E_out(E_in_array_edges, thickness)

        E_out_array_edges[np.isnan(E_out_array_edges)] = 0
        first_non_zero_energy = np.argmax(E_out_array_edges > 0.005) # since the stopping power below 5 keV isn't in the files

        dE_in_array = np.diff(E_in_array_edges)
        dE_out_array = np.diff(E_out_array_edges)
        dE_out_array[:first_non_zero_energy] = 1000000000

        yields_out_array = yields_in_array*dE_in_array/dE_out_array
        #yields_out_array[yields_out_array == np.infty] = 0
        #yields_out_array[first_non_zero_energy-1] = 0
        #yields_out_array[first_non_zero_energy] = 0
        #yields_out_array[first_non_zero_energy+1] = 0
        E_out_array = get_centered_from_edge(E_out_array_edges)

        return E_out_array, yields_out_array 

        
    def E_in_spectrum(self, E_out_array, yields_out_array, thickness):
        """
        Finds the spectrum before the particles passed through the material.

        Parameters
        ----------
        E_out_array : array of floats (numpy array)
            The energy array for the spectrum. MeV units.

        yields_out_array : array of floats (numpy array)
            The array of yields per MeV of the in spectrum.

        Returns
        -------
        E_in_array : array of floats (numpy array)
            The energy of the in spectrum. MeV units.

        yields_in_array : array of floats (numpy array)
            The array of yields per MeV of the in spectrum .

        """

        E_out_array_edges = get_edges_from_centered(E_out_array)
        E_in_array_edges = self.E_in(E_out_array_edges, thickness)

        dE_out_array = np.diff(E_out_array_edges)
        dE_in_array = np.diff(E_in_array_edges)

        yields_in_array = yields_out_array*dE_out_array/dE_in_array
        E_in_array = get_centered_from_edge(E_in_array_edges)

        return E_in_array, yields_in_array 
