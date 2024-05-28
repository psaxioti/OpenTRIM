# Iradina++

A C++ BCA Monte-Carlo code for performing SRIM-like simulations of ion
transport in materials.

The emphasis is on calculation of target damage.

## Usage

Iradina++ provides a command line program which can be invoked by 

```
> iradina++ [options] [-f config.json]
```
The program accepts a JSON-formatted configuration input either
directly from a file (with the `-f` option) or from stdin.

To see all available options run `iradina++ -h`, which prints
```
Monte-Carlo ion trasport simulation
Usage:
  iradina++ [OPTION...]

  -n arg          Number of histories to run (overrides config input)
  -j arg          Number of threads (overrides config input)
  -s, --seed arg  random generator seed (overrides config input)
  -f arg          JSON config file
  -v, --version   Display version information
  -h, --help      Display short help message
```

The program first checks and validates the configuration input. 

It then runs the simulation and saves the results into HDF5 files.

Additionally, the following libraries are provided, which can be used independently:
- `libiradinapp` contains all the Monte-Carlo ion transport C++ code as a shared library. It can be used by external applications
- `libiondedx` is a shared library containing tables of electronic stopping data
- The [`xs.h`](src/xs.h) C++ header file can be used standalone for screened Coulomb scattering calculations.
- `libcorteo4bit` & `libcorteo6bit` are shared libraries containing scattering tables of the ZBL potential. They can be used independently together with `xs.h` for scattering simulations.

## Input configuration

The JSON configuration input has the following self-explanatory structure:

```json
{
    "Simulation": {
        "simulation_type": "FullCascade",       // FullCascade|IonsOnly|CascadesOnly
        "screening_type": "ZBL",       // None|LenzJensen|KrC|Moliere|ZBL
        "scattering_calculation": "Corteo4bitTable", // Corteo4bitTable|Corteo6bitTable|ZBL_MAGICK|GCQuad
        "flight_path_type": "AtomicSpacing",    // AtomicSpacing|Constant|MendenhallWeller
        "straggling_model": "YangStraggling",   // NoStraggling|BohrStraggling|ChuStraggling|YangStraggling
        "nrt_calculation": "NRT_element",       // NRT_element|NRT_average
        "flight_path_const": 0.1,   // [nm], used if "flight_path_type"=Constant        
        "min_energy": 1.0,           // [eV], transport cutoff energy.
        "min_recoil_energy" : 1.0    // [eV], recoil energy cutoff, MendenhallWeller flight path selection   
    },
    "IonBeam": {
        "ion_distribution": "SurfaceCentered", // SurfaceRandom|SurfaceCentered|FixedPos|VolumeCentered|VolumeRandom
        "ionZ": 1,          // ion atomic number
        "ionM": 1.00784,    // ion mass
        "ionE0": 1e+06,     // ion energy [eV] 
        "dir": [1.0,0.0,0.0], // direction vector
        "pos": [0.0,0.0,0.0]  // [nm], used if "ion_distribution"=Fixed
    },
    "Target": {
        "materials": { 
            "Fe Layer": { // example material definition
                "density": 7.8658, // [g/cm^3]
                "Z": [26],      // atomic numbers
                "M": [55.845],  // atomic masses
                "X": [1.0],     // atomic fractions 
                "Ed": [40.0],   // displacement thresholds [eV]
                "El": [3.0],    // Lattice binding energies [eV] 
                "Es": [3.0],    // Surface binding energies [eV]
                "Er": [40.0]    // Replacement threshold energies [eV]
            } // more materials can be added here
        },
        "regions": { // rectangular regions filled with a material
            "R1": { // example region definition
                "material_id": "Fe Layer", // material must be defined above
                "min": [0.0,0.0,0.0],      // [nm] lower left corner 
                "max": [10000.0,10000.0,10000.0] // [nm] upper right corner
            }
        },
        "cell_count": [100,1,1], // # of cells in x,y,z-axis  
        "cell_size": [100.0,10000.0,10000.0], // cell width [nm] in x,y,z
        "periodic_bc": [0,1,1] // periodic boundary conditions in x,y,z
    },
    "Output": {
        "title": "Ion Simulation", // title written in output file
        "OutputFileBaseName": "iradina++", // base name for output files
        "store_transmitted_ions": 1, // store all ions that exit the simulation
        "store_pka": 1, // store all PKA events
        "store_dedx": 1 // store the electronic energy loss tables
    },
    "Driver": {
        "max_no_ions": 20000, // ion histories to run
        "threads": 8, // threads to use
        "seeds": [] // array of seeds, one per thread. If empty, random seeds are used
    }
}
```

> Comments in JSON are accepted.

The easiest way to get started is to copy/paste the above configuration make changes and save to a new input file.

Most of the options shown here are default values and can be omitted.
The `target/materials` and `target/regions` must always be given.

## Output files

The simulation produces a main file with the standard tallies and a number of
 optional files. 
 
 All output files will be named according to the config option `Output/OutputFileBaseName`.

 A brief description is given here. More detailed information can be found within the files. 

- `[basename].h5` is the main output file
- `[basename].pka.h5` [optional] contains a table with all PKA events
- `[basename].exit.h5` [optional] contains a table with all ions that left the simulation volume

The main output file has the following structure with the data organized in folders.

**Main output HDF5 file structure**
- `/` : Root folder
  - `Title` : the title given in the config
  - `Nh` : # of histories
  - `Variable_List` : Detailed list of variables in the file
  - `config_json` : JSON formatted configuration
  - `eels/` : electron energy loss tables folder
    - `dEdx_erg` : array of energy values [eV]
    - `dEdx` : stopping power [eV] table [energy x atoms x materials]    
    - `dEstrag` : straggling [eV/nm^(1/2)] table [energy x atoms x materials]
    - `maxImpactPar`
    - `max_fp`
  - `grid/` : geometric grid folder
    - `X` : x-axis grid points [nm]
    - `Y` : y-axis grid points [nm]
    - `Z` : z-axis grid points [nm]
    - `cell_xyz` : cell center positions [cells x 3]
  - `tally/`  : tallies folder
    - `totals/` : total defects/ions in the simulation volume
      - `Ndisp` : displacements
      - `Nimpl` : implantations / interstitals
      - `Nlost` : total # of ions that left the simulation
      - `Npkas` : PKAs
      - `Nrepl` : replacements
      - `Nvac`  : vacancies
    - `damage/` : damage quantities
      - `Tdam` : damage energy per pka [eV] [cells x atoms]
      - `Tdam_LSS` : LSS damage energy per pka [eV] [cells x atoms]
      - `Vnrt` : NRT vacancies from Tdam [cells x atoms]
      - `Vnrt_LSS` : NRT vacancies from Tdam_LSS [cells x atoms]
    - `defects/`
      - `Implantations` : Implanted ions and interstitials [cells x atoms]
      - `Replacements` : Replacements [cells x atoms]
      - `Vacancies` : Vacancies [cells x atoms]
      - `PKAs` : PKAs [cells x atoms] 
      - `Lost` : ions exiting the simulation [cells x atoms]
    - `energy_deposition/` : [eV]
      - `Ionization` : Ion energy deposited to ionization [cells x atoms]
      - `Phonons` : Ion energy deposited to phonons [cells x atoms]
      - `PKAs` : PKA recoil energy [cells x atoms]
      - `Lost`  : Energy lost due to ion exit [cells x atoms]
    - `ion_stat/` : Ion statistics
      - `flight_path` : ion path [nm] [cells x atoms]
      - `collisions` : # of ion collisions [cells x atoms]

To reach a variable in the file use the complete path, e.g. `/tally/damage/Tdam`.

The tallies give the mean values over all histories.
For each tally variable there is an additional entry corresponding to the standard deviation of the mean. The name of this entry is the same as the variable plus `_std` at the end, e.g.   `/tally/damage/Tdam_std`.

Specifically, if $x_i$ is the contribution to quantity $x$ from the $i$-th ion history, then the mean and std deviation given in the output file are calculated as:
$$
\bar{x} = \frac{1}{N_h} \sum_i {x_i}
$$
$$
\sigma_{\bar{x}}^2 = \frac{1}{N_h(N_h-1)} \sum_i { (x_i - \bar{x})^2 }
$$

## Code Documentation

Doxygen documentation can be found here: https://fusion.ipta.demokritos.gr/iradina++/

## Build & Install - Linux

The project can be built and installed with `cmake`.

The `Eigen` and `HDF5` C++ libraries are needed for building. Install them with
```bash
# Ubuntu 22.04 / DEB
sudo apt install libeigen3-dev libhdf5-dev libhdf5-103-1 libhdf5-cpp-103-1

# RHEL 9 / RPM
sudo dnf install eigen3-devel.noarch hdf5.x86_64 hdf5-devel.x86_64
```  

The HDF5 runtime libraries are needed for running the program.

Basic build receipe (run from project directory)

```
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make
make install
```
The default install location is `$HOME/.local`, thus `sudo` is not required.
Override this by setting the option `-DCMAKE_INSTALL_PREFIX="/your/install/location"` when calling `cmake`. 

Installed files are listed in `build/install_manifest.txt`, which is created in the install step.

## Testing

Some benchmark tests for comparison to other codes are given in folder `test/`.

The file [`test/README.md`](test/README.md) gives a short description of the benchmarks.

Folders [`test/iradina++/b1`](test/iradina++/b1) to `b7` have config files for running the benchmarks with iradina++.

The file [`test/octave/plot_benchmark.m`](test/octave/plot_benchmark.m) is an OCTAVE script which can be used for plotting results from benchmarks.

## Credits

Iradina++ draws heavily on the following similar open-source projects:

- The program [iradina](https://sourceforge.net/projects/iradina/) written by Ch. Borschel & C. Ronning and extended by J.P. Crocombette & Ch. Van Wambeke.
- The program [Corteo](http://www.lps.umontreal.ca/%7Eschiette/index.php?n=Recherche.Corteo) by F. Schiettekatte.

Electronic energy loss data has been obtained from the [SRIM-2013](http://www.srim.org/) distribution by  [J.F. Ziegler](ziegler[at]srim.org) using the provided utility `SRmodule.exe`

The [Eigen](http://eigen.tuxfamily.org/) library by B. Jacob & G. Guennebaud is used for 3D vector math.

The [Xoshiro256+](https://prng.di.unimi.it/) algorithm by D. Blackman and S. Vigna is used for random number generation.

[JSON for Modern C++](https://github.com/nlohmann/json) by N. Lohmann is used for encoding/decoding program options to/from json.

[cxxopts](https://github.com/jarro2783/cxxopts) by [jarro2783](https://github.com/jarro2783) is used for handling cli options.

The [HDF5 library](https://github.com/HDFGroup/hdf5) is used for saving
results to HDF5 files.