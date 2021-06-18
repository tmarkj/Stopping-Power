#!/usr/bin/python3

import stopping_power
import numpy as np

H_in_Ta = stopping_power.Stopping_power('H', 'Ta')
H_in_Ni = stopping_power.Stopping_power('H', 'Ni')

E_in = 3.0 #MeV

E_out_Ta = H_in_Ta.E_out(E_in, 10)

E_out_Ni = H_in_Ni.E_out(E_out_Ta, 67)

print(E_out_Ta)
print(E_out_Ni)

## Testing array inputs

E_in  = np.linspace(.4, 4, 10)

E_out_array_Ta = H_in_Ta.E_out(E_in, 10)
print(E_out_array_Ta)

print(H_in_Ta.range(20.0))
print("3.3", H_in_Ta.range(3.3))
print("3.314", H_in_Ta.range(3.314))
print("3.35", H_in_Ta.range(3.35))
print("E_in 3.5 thickness 10: ", H_in_Ta.E_out(3.5, 10))

H_in_Al = stopping_power.Stopping_power('H', 'Al')
print("Aluminum 3.5", H_in_Al.range(3.5))

E_in_array = [2.8,2.9,3.0,3.1,3.2]
yield_in_array = [.25,.75,1.5,.75,.25]

print(H_in_Al.E_out_spectrum(E_in_array, yield_in_array, 75.0))
