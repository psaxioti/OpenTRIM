# Iradina++

A C++ BCA Monte-Carlo code for performing SRIM-like simulations of ion
transport in materials.

The emphasis is on calculation of target damage.

## Usage

Iradina++ is a command line program and can be invoked by 

```
> iradina++ [options] [ < [config_file.json]]
```
The program accepts a JSON-formatted configuration either
directly from stdin or from a file `[config_file.json]` by redirection
as shown above.

Currently the only option that can be given is 
```
  -h : print a short help message and exit.
```

The program first checks and validates the input. If something is wrong
it stops and indicates the error.

It then runs the simulation and saves the results into HDF5 files.

## Input configuration

The JSON configuration input has the following self-explanatory structure:

```json
{
    "Simulation": {
        "simulation_type": "FullCascade",       // FullCascade|IonsOnly|CascadesOnly
        "nrt_calculation": "NRT_element",       // NRT_element|NRT_average
        "scattering_calculation": "Corteo4bit", // Corteo4bit|Corteo6bit|ZBL_MAGICK
        "flight_path_type": "AtomicSpacing",    // AtomicSpacing|Constant|MendenhallWeller
        "straggling_model": "YangStraggling",   // NoStraggling|BohrStraggling|ChuStraggling|YangStraggling
        "flight_path_const": 0.1,   // [nm], used if "flight_path_type"=Constant
        "min_energy": 1.0           // [eV], min cutoff ion energy 
    },
    "IonBeam": {
        "ion_distribution": "SurfaceCentered", // SurfaceCentered|SurfaceCentered|FixedPos|VolumeCentered|VolumeRandom
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
                "M": [55.845],  // atomic mass
                "X": [1.0],     // atomic concentration (normalized by the program)
                "Ed": [40.0],   // displacement threshold [eV]
                "El": [3.0],    // Lattice binding [eV] 
                "Es": [3.0],    // Surface binding [eV]
                "Er": [40.0]    // Replacement threshold [eV]
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

> Typically comments are not allowed in JSON but for iradina++ it is ok. 

Copy/paste and edit the above to create a new input file.
Most of the options shown here are default values and can be omitted.
The target materials and regions must always be given.

## Output files

The simulation produces a main file with the standard tally output and a number of
 optional files. All output files will be named according to the config option Output/OutputFileBaseName.

 A brief description is given here. Detail information can be found within the files 

- `[basename].h5` is the main output file that contains:
  - The input configuration
  - Cell coordinates
  - Tally data for generated vacancies, implanted ions, interstitials, replacements, lost ions 
  - Tally data for energy deposition in electronic ionization, the lattice and lost due to ions leaving the simulation
  - Tally data for PKA damage energy and NRT estimated vacancy generation
  - [optional] Tables of electronic energy loss and straggling
- `[basename].pka.h5` [optional] contains a table with all PKA event
- `[basename].exit.h5` [optional] contains a table with all ions that left the simulation volume

## Documentation

Doxygen documentation can be found here: https://fusion.ipta.demokritos.gr/iradina++/

## Credits

Iradina++ draws heavily on the following similar open-source projects:

- The program [iradina](https://sourceforge.net/projects/iradina/) written by Ch. Borschel & C. Ronning and extended by J.P. Crocombette & Ch. Van Wambeke.
- The program [Corteo](http://www.lps.umontreal.ca/%7Eschiette/index.php?n=Recherche.Corteo) by F. Schiettekatte.

Electronic energy loss data has been obtained from the [SRIM-2013](http://www.srim.org/) distribution by  [J.F. Ziegler](ziegler[at]srim.org) using the provided utility `SRmodule.exe`

The [Xoshiro128+](https://prng.di.unimi.it/) algorithm by D. Blackman and [S. Vigna](vigna@acm.org) is used for random number generation.

[JSON for Modern C++](https://github.com/nlohmann/json) by N. Lohmann is used for encoding/decoding program options to/from json.

The [HDF5 library](https://github.com/HDFGroup/hdf5) is used for saving
results to HDF5 files.