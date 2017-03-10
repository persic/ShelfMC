# ShelfMC
Monte Carlo for simulating the sensitivity of UHE neutrino detectors in ice.

Important note: input_reference.txt includes an absolute file location which is specific to Chris's installation on the UCI HPC cluster, and must be modified.  Similarly, the scripts to generate qsub files are also specific to HPC, and must be customized for your cluster.

Station types are defined as follows:
0. Isotropic antennas (antenna type 0) are spaced equally around the station center. Variable number of antennas.
1. This is the pre-2017 shelfMC station.  LPDA's (type 1) are placed equally around the station center, polarized in the tangential direction, facing down. Variable number of antennas.
2. LPDA's based on the response from WIPL-D simulations(type 2) are placed equally around the station center, polarized in the tangential direction, facing down. Variable number of antennas.
3. 8 LPDA's based on the response from WIPL-D simulations(type 2) are placed in pairs, facing out from the center.  All antennas are H-pol.
4. 8 LPDA's based on the response from WIPL-D simulations(type 2) are placed in	pairs, facing out from the center.  All	antennas are V-pol.
 

Antenna types are defined as follows
0. Totally isotropic antenna.  No direction, polarization, or frequency dependence on gain.
1. Pre-2017 LPDA model.  Gaussian directional response.
2. 100MHz Create LPDA model from WIPL-D simulations. 

