stopping_power.so : stopping_power.c
	gcc -O2 -shared -o stopping_power.so -fPIC stopping_power.c
