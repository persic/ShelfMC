# ShelfMC

Monte Carlo for simulating the sensitivity of UHE neutrino detectors in ice.

Getting Started
---------------
### Pre-Reqs

  1. ShelfMC compiles using ROOT libraries, and saves output using the ROOT tree data structure, so you'll need to get yourself some ROOT.  ROOT version 5.34 has been proven compatible.  Root 6 has not been tested.  If you need to install ROOT, follow the directions here https://root.cern.ch/building-root.

  2. ShelfMC is written in C++, but some of the supplementary scripts are python based, and require Numpy and/or Matplotlib libraries.  See https://www.scipy.org/ for more info.

  3. Certain antenna types require access to certain files which define the antenna response.  These files are too big for GitHub, so you need to contact Chris or Anna, and they'll get them to you.  If you don't know who Chris or Anna are, ask the person who told you to run ShelfMC in the first place.  These files should be placed in the main ShelfMC directory.

### Building ShelfMC

To build ShelfMC, navigate to your installation directory and run

```bash
  make
```
Hopefully this works, and you now have an executable file named ```shelfmc_stripped.exe```

### Running ShelfMC

  1. First create a working directory where ShelfMC will save its output.  Before running ShelfMC, this working directory must contain a file called ```input.xml``` which contains the program's input parameters.  This file should be based off the ```input_reference.xml``` file which lives in the main ShelfMC directory, which you can reference for definitions of the various input parameters.

  2. A few important notes about the working environment for ShelfMC.
    * Relative file paths are based on the directory ShelfMC is run from, NOT the working directory that contains the input.xml file.
    * The antenna model files mentioned in "Pre-Reqs" are defined with relative file paths in ```shelfmc_stripped.cc```, so make sure they are in the directory where you run ShelfMC.



  3. Once your running environment and working directory are set up, you can execute ShelfMC by running the command

  ```bash
  ./shelfmc_stripped.exe [path/to/workingDirectory/] [Unique_Tag_For_Your_Output_Files]
  ```
  You may wish to direct the output into a file for later reference, or better yet, write some code to execute multiple instances of ShelfMC on your friendly neighborhood super-computer.

Tips & References
-----------------

### Antenna Types

Currently, the various antenna types that you can use in ShelfMC are defined in the source code (see the GetHeff function in ```functions.cc```), so if you want to add a new antenna type, you'll need to do it there.  The current antenna types supported in ShelfMC are...

  0. Totally isotropic antenna.  No direction, polarization, or frequency dependence on gain.
  1. Pre-2017 Theoretical LPDA model.  Gaussian directional response.
  2. 100MHz Create LPDA model from WIPL-D simulations.
  3. Simulated ARA Bicone dipole antenna

### Station Geometry

The position, orientation, and type of the antennas that make up a station (for triggering purposes) are defined in an .xml document.  This is easily customizable, and examples can be found in the ```StnGeoFiles/``` directory.
