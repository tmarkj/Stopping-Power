#!/usr/bin/python3

import stopping_power

H_in_Ta = stopping_power.Stopping_power('H', 'Ta')
H_in_Ni = stopping_power.Stopping_power('H', 'Ni')

E_in = 3.0 #MeV

E_out_Ta = H_in_Ta.E_out(E_in, 10)

E_out_Ni = H_in_Ni.E_out(E_out_Ta, 67)

print(E_out_Ta)
print(E_out_Ni)
