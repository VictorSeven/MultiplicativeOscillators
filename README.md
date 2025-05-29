# MultiplicativeOscillators

Code repository for [Buendía PRL 2025](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.134.197201). If the code is used, please consider also citing the publication.

## Organization of the code

There are three main folders. The user interested to run the code should focus on the Julia version always. However, I leave the other ones as they contain code I used for the original paper (replicability).

- `cpp` contains the code for the microscopic Kuramoto model in C++. I did many simulations of the microscopic model using this code, but it got displaced by Julia later on. You can find more information on this code at the end of this document.
- `python` contain launchers for the C++ code. 
- `julia` contains the bulk of the code for microscopic and mesoscopic model simulation, as well as the paper figures. 

## Julia code

### Core computations

This is the main code used for the above publication. It contains many useful functions, particularly regarding the mesoscopic equations, that can be re-used. The `core` folder contains several files with the functions needed to simulate the system. Each file is a different version of the model and comes with a different module. The functions are called the same in each one and return the same information so it is possible to change the integrated model by calling the correct module. 

The files are the following:

- `amplitude_additive.jl`. Integrates all the equations subject to independent additive noise and phase ansatz. It includes finite-size corrections only to the Kuramoto parameter, others are approximated to zero to ensure numerical stability. Commented in the file there is included also a multiplicative noise version which I also employed to test, but the results seem to be largely independent. The module is called `AdditiveNoise`.

- `amplitude_stochastic.jl`. Integrates the full system using the phase ansatz. It includes finite-size corrections only to the Kuramoto parameter, others are approximated to zero to ensure numerical stability. The module is called `AmplitudeEquations`.

- `reduced-kuramoto-cartesian.jl`. Integrates the full system using the `x,y` variables instead of amplitudes and phases. Behaves in a more stable way, taking into account all finite sizes without a problem, but it is slower -as it is double the times of equations. The module name is `CartesianEquations`.

- `amplitude_fixedpointnoise.jl`. Not used in the paper but included for convenience. It uses the linear noise approximation, which is very accurate close to an equilibrium point. Here, the correlated multiplicative noise is evaluated once before the simulation starts and it is considered constant throughtout. Uses the phase ansatz and finite size corrections to Kuramoto order parameter only. The module name is `FixedPointNoise`.

- `kuramoto.jl` contains the integration of the microscopic mean-field Kuramoto model. The module is called `Kuramoto`.

All of these files have the same functions for ease of use. The most important are

- `get_timeseries(N, w, q, s2, trelax, tf, filename; dt=0.01, nsample=100, nharms=1)` computes a timeseries for the given parameters. It lets the system relax for `trelax` and then integrates for `[0,tf]`. It stores the results in `filename`, writing each `nsample` iterations, and stores the specified number of harmonics. It writes a file with `t r1 r2 ... rnharms` at each row.

- `phase_diagram(N, w, q0, qf, nq, s2, trelax, tf, filename; dt=0.01, nsample=100)`. It computes the average Kuramoto parameter in the interval `[q0, qf]` using `nq` steps. Most parameters are identical to the previous function. Here `nsample` is the number of iterations left between measurement of the R,psi values. THis avoids correlations between them. It writes a file with three columns, `q, av_r, var_r` at each row.

Finally, the `theory_formulas.jl` file contains several theoretical formulas and models used to compare.

### Launchers

I used a computing cluster facility to perform the computations. The `launcher` folder contains `slurmscript` files that call Julia files that perform the computations. These Julia file setup the system and call the `core` functions mentioned earlier. They can serve as an example to reuse the code. 

### Figure replication

In order to replicate the figures of the papers, the data would be necessary. This data is heavy and not provided, but it can be generated from the user by using the rest of the code. Then, it suffices to just call `julia figureX.jl` to produce the PDF files seen in the paper (as well as some larger-sized figures I use for the talks). 


## C++ / Python code

This code was used for a long time to produce simulations for the microscopic system. The `kuramoto-mf.cpp` file works by selecting a mode at compile time. For timeseries, the number of harmonics to proudce is also set at compile time, like `g++ -o3 -DMODE=SERIES -ORDER=10`. Use `-DMODE=DIAGRAM` to generate phase diagrams. For the `DIAGRAM` mode, the program takes the arguments `N,w,s,q0,qf,nq,trelax,tf,filename,seed`, similar to the Julia function. In the `SERIES` mode, it needs `N,w,s,q,trelax,tf,filename`. The launchers for the C++ code are written in Python. There is no integrator for the reduced model. In general, for new code the Julia version should be preferred, but the code is left here for archival and replicability purposes.


