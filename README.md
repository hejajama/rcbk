Solve the BK equation in coordinate space with running coupling

Questions and comments: Heikki Mäntysaari <heikki.mantysaari@jyu.fi>

License: CC-BY

Reference: T. Lappi, H. Mäntysaari, Phys.Rev.D 88 (2013) 114020, arXiv:1309.6963 

Parallerized using OpenMP

Compile (requires GSL and CMake)
```bash
mkdir build
cd build
cmake ..
make
```

The datafiles containing the solution (evolved dipole amplitude) can be read using the code [rcbkdipole](https://github.com/hejajama/rcbkdipole)

## Examples
Solve up to $y=10$, save data to `datafile`, use the Balitsky running coupling prescription, set $C^2=10$, initial condition at $x_0=0.01$ is the MV model with $Q_{s,0}^2=1$ GeV$^2$, $\gamma=1.2$, $e_c=2$.
```bash
./build/bin/rcbk -ic MV 1 1.2 0.01 2 -alphas_scaling 10 -rc BALITSKY -output datafile -maxy 10 -minr 1e-6 -fast 2>/dev/null
```
Notes
* Without the `-fast` flag the code is very slow, and the effect at the level of the dipole amplitude is negligible, except maybe if you want to compute the Fourier transform
* Lower limit of the dipole size grid is here set to $r_\mathrm{min}=10^{-6}$ $\mathrm{GeV}^{-1}$. Especially at small $r$ the required numerical precision may not be reached
* and the code may print a lot of warnings, thus in the example above these warnings are directed to `/dev/null` (=hidden)

Same example, but read the initial condition at $x_0=0.01$ from a datafile where the structure is 
```
r1 N(r1)
r2 N(r2)
...
rN N(rN)
```
(here `ri` is in GeV $^{-1}$ )

```bash
./build/bin/rcbk -ic FILE initial_condition_file 0.01 -alphas_scaling 10 -rc BALITSKY -output datafile -maxy 10 -fast 2>/dev/null
```

Read the `main.cpp` to see all other parameters
