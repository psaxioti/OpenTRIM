# User Documentation {#userdoc}

- @ref install
- @ref cliapp
- @ref json_config
- @ref out_file

## Installation  {#install}

### Linux

The project can be built and installed with `cmake`.

The `Eigen` and `HDF5` libraries are needed for building. Install them with
```bash
# Ubuntu 22.04 / DEB
sudo apt install libeigen3-dev libhdf5-dev libhdf5-103-1

# RHEL 9 / RPM
sudo dnf install eigen3-devel.noarch hdf5.x86_64 hdf5-devel.x86_64
```  

The HDF5 runtime libraries are needed for running the program.

Basic build recipe (run from project directory)

```
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make
make install
```
The default install location is `$HOME/.local`, thus `sudo` is not required.
Override this by setting the option `-DCMAKE_INSTALL_PREFIX="/your/install/location"` when calling `cmake`. 

### Windows

Download the latest binary distribution release as a zip file and extract to some location.

To run the program from the windows command line, add the bin/ subfolder to the system or user path.

## Using the command line application {#cliapp}

Iradina++ provides a command line program which can be invoked by 

    iradina++ [options] [-f config.json]

The program accepts a @ref json_config "JSON-formatted configuration input" either
directly from a file (with the `-f` option) or from stdin.

To see all available options run `iradina++ -h`, which prints

```
Monte-Carlo ion trasport simulation
Usage:
  iradina++ [OPTION...]

 -n arg            Number of histories to run (overrides config input)
 -j arg            Number of threads (overrides config input)
 -s, --seed arg    random generator seed (overrides config input)
 -o, --output arg  output file base name (overrides config input)
 -f arg            JSON config file
 -t, --template    pring a template JSON config to stdout
 -v, --version     Display version information
 -h, --help        Display short help message
```

The program first checks and validates the configuration input. 

It then runs the simulation and saves the results into a HDF5 archive.

## The JSON configuration file {#json_config}

All configuration parameters for a simulation can be coded in a JSON formatted string. 

This can be loaded directly from a file to the \ref cliapp "iradina++ cli program" with the following command:

    iradina++ -f config.json

where `config.json` is a file containing the JSON configuration.

In a C++ program one can use the \ref mcdriver::options class which offers
a parseJSON() function to parse and validate the options.

The JSON-formatted options string has the following self-explanatory structure shown bellow.
This example has all default options and an example target definition.

```javascript
{
    "Simulation": {
        "simulation_type": "FullCascade",       // FullCascade|IonsOnly|CascadesOnly
        "screening_type": "ZBL",       // None|LenzJensen|KrC|Moliere|ZBL
        "scattering_calculation": "Corteo4bitTable", // Corteo4bitTable|ZBL_MAGICK
        "flight_path_type": "AtomicSpacing",    // AtomicSpacing|Constant|MendenhallWeller|IPP
        "straggling_model": "YangStraggling",   // NoStraggling|BohrStraggling|ChuStraggling|YangStraggling
        "nrt_calculation": "NRT_element",       // NRT_element|NRT_average
        "flight_path_const": 0.1,   // [nm], used if "flight_path_type"=Constant        
        "min_energy": 1.0,           // [eV], transport cutoff energy.
        "min_recoil_energy" : 1.0,    // [eV], recoil energy cutoff, flight_path_type = MendenhallWeller or IPP 
        "allow_sub_ml_scattering": false, // allow flight path smaller than 1ML, flight_path_type = IPP
        "max_mfp": 3.4028235e+38, // upper bound for ion mfp, flight_path_type = IPP
        "max_rel_eloss": 0.05 // upper bound of relative electronic energy loss, flight_path_type = MendenhallWeller or IPP   
    },
    "IonBeam": {
        "ion_distribution": "SurfaceCentered", // SurfaceRandom|SurfaceCentered|FixedPos|VolumeCentered|VolumeRandom
        "ionZ": 1,          // ion atomic number
        "ionM": 1.00784,    // ion mass
        "ionE0": 1e+06,     // ion energy [eV] 
        "dir": [1.0,0.0,0.0], // direction vector
        "pos": [0.0,0.0,0.0]  // [nm], used if "ion_distribution"=FixedPos
    },
    "Target": {
        "materials": { 
            "Bulk Iron": { // example material definition
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
                "material_id": "Bulk Iron", // material must be defined above
                "min": [0.0,0.0,0.0],      // [nm] lower left corner 
                "max": [10000.0,10000.0,10000.0] // [nm] upper right corner
            }
        },
        "cell_count": [100,1,1], // # of cells in x,y,z-axis  
        "cell_size": [100.0,10000.0,10000.0], // cell width [nm] in x,y,z
        "periodic_bc": [0,1,1] // periodic boundary conditions in x,y,z
    },
    "Output": {
        "title": "Ion Simulation", // short title for the simulation 
        "OutputFileBaseName": "out", // base name of output file (extension .h5)
        "store_transmitted_ions": 1, // store all ions that exit the simulation
        "store_pka": 1, // store all PKA events
        "store_dedx": 1 // store the electronic energy loss tables
    },
    "Driver": {
        "max_no_ions": 20000, // ion histories to run
        "threads": 8, // threads to use
        "seeds": [] // seed to use for random number generation. If empty, a random seed is generated
    }
}
```

@note Iradina++ accepts comments in the JSON config string

The easiest way to get started is to get a template with all default options by running

    iradina++ -t > template.json

and make changes to the new template file.

Most of the options shown above are default values and can be omitted.

Only the `target/materials` and `target/regions` must always be given.

## The HDF5 output archive {#out_file}

The simulation produces a single HDF5 output file which contains all results from tallies and events along with the input conditions, run timing and other information. 
 
A brief description of the file contents is given here. More detailed information can be found within the file itself. 

Dimensions of tables depend on:
- \f$N_{at}\f$ : # of atoms, i.e., all atoms in the target plus the projectile. Note that the same atom species belonging to different materials is counted as different atom. 
- \f$N_{mat}\f$ : # of target materials
- \f$N_x, N_y, N_z\f$ : # of grid points along the 3 axes
- \f$N_c=(N_x-1)(N_y-1)(N_z-1)\f$ : # of target cells
- \f$N_e\f$ : # of energy points for energy loss tables
- \f$N_{ev}\f$ : # of events

**HDF5 output archive structure** 

- `/` Root folder
  - `config_json`: [1x1], JSON formatted simulation options (string)
  - `title`: [1x1], user supplied simulation title (string)
  - `variable_list`: [1x1] list of variables in file (string)
  - `atom/`
    - `Ed`: [\f$N_{at}\times 1\f$], displacement energy [eV]
    - `El`: [\f$N_{at}\times 1\f$], lattice binding energy [eV]
    - `Er`: [\f$N_{at}\times 1\f$], replacement energy [eV]
    - `Es`: [\f$N_{at}\times 1\f$], surface binding energy [eV]
    - `M`: [\f$N_{at}\times 1\f$], atomic mass
    - `Z`: [\f$N_{at}\times 1\f$], atomic number
    - `label`: [\f$N_{at}\times 1\f$], labels = [Atom (Chemical name)] in [Material]
    - `name`: [\f$N_{at}\times 1\f$], Chemical names
  - `eels/`
    - `dEdx`: [\f$N_{at}\times N_{mat}\times N_e\f$], dEdx values [eV/nm]
    - `dEdx_erg`: [\f$N_e \times 1\f$], dEdx table energy values [eV]
    - `dEstrag`: [\f$N_{at}\times N_{mat}\times N_e\f$], straggling values [eV/nm^(1/2)]
    - `ipmax`: [\f$N_{at}\times N_{mat}\times N_e\f$], max impact parameter [nm]
    - `mfp`: [\f$N_{at}\times N_{mat}\times N_e\f$], ion mean free [nm]
  - `exit_events/`
    - `column_descriptions`: [10x1], event data column descriptions
    - `column_names`: [10x1], event data column names
    - `event_data`: [\f$N_{ev}\times 10\f$]
  - `grid/`
    - `X`: [\f$N_x\times 1\f$], x-axis grid
    - `Y`: [\f$N_y\times 1\f$], y-axis grid
    - `Z`: [\f$N_z\times 1\f$], z-axis grid
    - `cell_xyz`: [\f$3\times N_c\f$], cell center \f$(x,y,z)\f$ coordinates
  - `material/`
    - `atomic_density`: [\f$N_m\times 1\f$], atomic density [at/nm^3]
    - `atomic_radius`: [\f$N_m\times 1\f$], atomic radius [nm]
    - `mass_density`: [\f$N_m\times 1\f$], mass density [g/cm^3]
    - `name`: [\f$N_m\times 1\f$], names of materials
  - `pka_events/`
    - `column_descriptions`: [8x1], event data column descriptions
    - `column_names`: [8x1], event data column names
    - `event_data`: [\f$N_{ev}\times 8\f$]
  - `run_stat/`
    - `Nh`: [1x1], # of histories
    - `end_time`: [1x1], finish time/date (string)
    - `ips`: [1x1], ion histories per second
    - `process_cpu_time`: [1x1], total cpu time [s]
    - `start_time`: [1x1], start time/date (string)
  - `tally/`  The score tables from the \ref tally object
    - `damage/`
      - `Tdam`: [\f$N_{at}\times N_c\f$], Damage energy [eV]
      - `Tdam_LSS`: [\f$N_{at}\times N_c\f$], Damage energy estimated by the LSS approximation [eV]
      - `Vnrt`: [\f$N_{at}\times N_c\f$], Vacancies per the NRT model using Tdam
      - `Vnrt_LSS`: [\f$N_{at}\times N_c\f$], Vacancies per the NRT model using Tdam_LSS
    - `defects/`
      - `Implantations`: [\f$N_{at}\times N_c\f$], Implantations
      - `Lost`: [\f$N_{at}\times N_c\f$], Ions that exit the simulation volume
      - `PKAs`: [\f$N_{at}\times N_c\f$], PKAs
      - `Replacements`: [\f$N_{at}\times N_c\f$], Replacements
      - `Vacancies`: [\f$N_{at}\times N_c\f$], Vacancies
    - `energy_deposition/`
      - `Ionization`: [\f$N_{at}\times N_c\f$0], Energy deposited to electron ionization [eV]
      - `Lost`: [\f$N_{at}\times N_c\f$], Energy lost due to ions exiting the simulation [eV]
      - `PKA`: [\f$N_{at}\times N_c\f$], PKA recoil energy [eV]
      - `Phonons`: [\f$N_{at}\times N_c\f$], Energy deposited to the lattice [eV], Phonons = Stored+Lattice
      - `Stored`: [\f$N_{at}\times N_c\f$], Energy stored in lattice defects [eV]
      - `Lattice`: [\f$N_{at}\times N_c\f$], Thermal energy deposited to the lattice [eV]
    - `ion_stat/`
      - `collisions`: [\f$N_{at}\times N_c\f$], ion collisions
      - `flight_path`: [\f$N_{at}\times N_c\f$], flight path [nm]
    - `totals/`
      - `data`: \f$[19\times 1]\f$ for each tally, sum of scores from all atoms and all cells
      - `column_names`: name of each column in `data`	
  - `version_info/`
    - `build_system`: [1x1], build operating system
    - `build_time`: [1x1], build timestamp
    - `compiler`: [1x1], compiler id
    - `compiler_version`: [1x1], compiler version
    - `version`: [1x1], iradina++ version

To reach a variable in the archive use the complete path, e.g. `/tally/damage/Tdam`.

Each number in the tally tables is a mean value over all histories.
For each tally table the standard error of the mean (SEM) is also included. This is a separate table of equal dimension and with the same name plus the ending `_sem`. E.g., for the table `/tally/damage/Tdam` there is also `/tally/damage/Tdam_sem`.

The definition of SEM employed is as follows: if \f$x_i\f$ is the contribution to the table entry \f$x\f$ from the \f$i\f$-th ion history, then the mean, \f$\bar{x}\f$, and the SEM, \f$\f$\sigma_{\bar{x}}\f$\f$, listed in the output file are calculated as:
$$
\bar{x} = \frac{1}{N_h} \sum_i {x_i}
$$
$$
\sigma_{\bar{x}}^2 = \frac{1}{N_h(N_h-1)} \sum_i { (x_i - \bar{x})^2 }=
\frac{\bar{{x^2}} - \bar{x}^2}{N_h-1}
$$

\see tally

