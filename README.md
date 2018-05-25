# SeAN
## Self Absorption Numerical

## Table of Contents

1. [Introduction](#introduction)
2. [Installation](#installation)

    2.1 [Prerequisites](#prerequisites)

    2.2 [Compilation](#compilation)

    2.3 [Build options](#build_options)

3. [Usage](#usage)

    3.1 [Command line options](#command_line_options)

    3.2 [Input file: Overview](#input_file_overview)

    3.3 [Input file: Options](#input_file_options)

4. [License](#license)

5. [References](#references)

### 1. Introduction <a name="introduction"></a>

`SeAN` is a program to numerically calculate the 'probability' (or, more precisely, the cross section) for the resonant absorption of photons in atomic nuclei.
The probability for absorption depends on the type of atomic nucleus and the system of atoms in which the nucleus is located. A special focus of `SeAN` is on the computationally efficient and correct inclusion of these 'condensed matter effects'.
`SeAN` can be used to simulate the effect of resonant absorption on a beam of photons passing through matter. It derives its name from a particular type of photonuclear experiment that can be used to measure the lifetime of excited states of a nucleus.
This README file only describes the usage of `SeAN` for someone who is already familiar with the subject. For a short introduction into the theory of photoabsorption on nuclei and the influence of condensed matter effects, see the documentation in `doc/` and references therein.

### 2. Installation <a name="installation"></a>

#### 2.1 Prerequisites <a name="prerequisites"></a>

Essential requirements:

* C++11 or more recent
* [ROOT](https://root.cern.ch/) (plotting, interpolation, integration) [tested with version 6.10/06]
* [FFTW](http://www.fftw.org/) (FFT algorithm for fast convolutions) [tested with version 3.3.7]

Optional:

* [openmp](http://www.openmp.org/) (parallelization)
* latex, bibtex (building the documentation)
* python3, matplotlib, scipy, numpy (running one of the unit tests)

#### 2.2 Compilation <a name="compilation"></a>

Get the latest version of `SeAN` by cloning it from the repository:

```
$ git clone https://github.com/uga-uga/SeAN.git
```

Enter the newly created directory `sean/` and execute

```
$ cd sean/
$ make
```

This should create a `sean` executable which is ready to run. Test it by executing 

```
$ ./sean --help
Usage: sean [OPTION...] INPUTFILE
SeAN, Self-Absorption Numerical
[...]
```

To remove all files which are created during the compilation, including the executable, type

```
$ make clean
```

#### 2.3 Build options <a name="build_options"></a>

At the moment, build options can only be set by manually changing the `sean/Makefile` and executing `make` again in the same directory.

##### No openmp

If you don't have openmp installed in your system, remove the `-fopenmp` compiler flag. This will result in warning messages during the compilation, but the program will run as usual (only slower):

```
[...] warning: ignoring #pragma omp [...]
```

##### Optimization

By default, the program is built as release version with the optimization flag `-O3`. To build in debug mode, use the `-g` option.

##### Change C++ compiler

If you don't have the default compiler `g++` or want to use another one on purpose, change the value of the `CC` variable in the Makefile.

### 3. Usage <a name="usage"></a>

`SeAN` reads all of its physics input parameters from an input file `INPUTFILE`. However, some command line options exist to control the output and the verbosity. They will be described in section [Command line options](#command_line_options). After that follows an introduction to the input file structure in section [Input file: Overview](#input_file_overview). [Input file: Options](#input_file_options) lists all possible options in the input file.

The general command to run `SeAN` is

```
$ ./sean [OPTION...] INPUTFILE
```

All output will be saved to the directory `output/`. Several other directories with self-explanatory names exist from which input is read.

#### 3.1 Command line options <a name="command_line_options"></a>

`SeAN` accepts the following command line options:

 * `-e`, `--exact`: Instead of using the FFT approximation to calculate the broadening of the Breit-Wigner cross section, evaluate the integral exactly. The speed advantage of using FFT is enormous, but the user is advised to check with this option whether it is applicable.
 * `-o`, `--output=OUTPUTFILENAME`: Write the input from `INPUTFILE` and the results of the calculation (the amount of resonant scattering on each of the targets) to a file called OUTPUTFILENAME.
 * `-p`, `--plot`: Create plots of all calculated intermediate quantities that can be shown in a histogram or graph. Note that this may take a lot of time and create very large files depending on the values of `NBINS_E` and `NBINS_Z`.
 * `-r`, `--recoil`: Consider the recoil of the nucleus when the absorption is calculated. This is a little tricky, because the absorption line is very narrow in the sub-eV range, but depending on the nucleus the recoil may shift the absorption line by up to several keV. This has to be considered when setting the integration range. The resulting resonant scattering is actually independent of the recoil, therefore it is left to the user to include this effect if strictly realistic results are desired.
 * `-v`, `--verbosity=VERBOSITY`: Set the verbosity on the command line.
   * `VERBOSITY == 0`: Print nothing
   * `VERBOSITY == 1`: Print results
   * `VERBOSITY == 2`: (default) Print input and results 
 * `-w`, `--write`: Create text files of all calculated intermediate quantities that can be shown in a histogram or graph. Note that this may take a lot of time and create very large files depending on the values of `NBINS_ENERGY` and `NBINS_Z`. This option can be used to avoid calculating the same thing (for example, an absorption cross section) over and over again, because the output can be read as input by the next calculation.
 * `-h`, `--help`: Print help to command line

#### 3.2 Input file overview <a name="input_file_overview"></a>

The structure of the input file is very analog to a nuclear physics experiment. One input file can be considered a single run in an experiment in which the physical system and the incoming photon intensity distribution is constant. The beam is supposed to traverse an array of materials, so-called 'targets'. `SeAN` will calculate the necessary properties of all targets, for example the absorption cross section, first, and then propagate the beam through the array.
Consequently, the input file consists of some general options and then an arbitrary amount of targets. The file is read line by line, ignoring the comment lines that start with a `#`.
In the following, the necessary input for a minimal calculation will be described. The examples which are referred to here can be found in `examples/basic/`. They describe a realistic self-absorption experiment (consisting of two separate runs) using one of the strongest gamma-ray transitions of all nuclei and the first gamma ray in the nuclear chart: the second excited state of 6Li.

##### General options

Before targets can be defined, the input file MUST start with four lines containing information that concerns all targets:

 * `EMIN, EMAX` (double, double): Minimum and maximum energy for the calculation in eV
 * `INCIDENT_BEAM , {PARAMETERS}` (string, ...): The model for the intensity distribution of the primary photon beam plus the parameters that this model needs.
 * `NBINS_ENERGY` (int): The number of bins on the energy axis
 * `NBINS_Z` (int): The number of bins on the z-axis, which is the axis in real space along which the beam propagates through the targets

The example files show a sane choice for the parameters. The absorption line is completely covered along with its slowly decaying tail. A decent number of bins on the energy- and z-axis has been chosen for a short execution time, even when plotting and writing. The primary beam is constant over the whole energy range with a dimensionless intensity of 1.

```
# Minimum and maximum energy for the calculation in eV
3.561880e6, 3.563880e6
# Primary beam intensity distribution
const, 1.0
# Number of bins on the energy axis
1000
# Number of bins on the z-axis
1000
```

##### Targets

Each target has 11 properties that need to be set. Since there is no introductional keyword for the properties, they all need to be set explicitly and exactly in the order shown below:

 * `IDENTIFIER` (string): Identifier string of the target, that will be used in all output. Its most important function is to define the names of the output files, for example 'IDENTIFIER_crosssection_0.txt'. Each target must have a unique identifier, otherwise files will be overwritten.

The following options are properties of the nuclear resonances. An arbitrary number of resonances can be defined (indicated by curly brackets), separated by commas. Obviously, the first resonance energy in the line of resonance energies belongs to the first spin quantum number on the line of spins and so on...

 * `{RESONANCE_ENRGIES}` ({double}): Resonance energies in eV
 * `{GAMMA_0}` ({double}): Transition widths to the ground state in eV
 * `{GAMMA}` ({double}): Total transition widths in eV
 * `J_0` (double): Angular momentum quantum number of the ground state in units of the reduced Planck constant
 * `{J_0}` ({double}): Angular momentum quantum numbers of the excited states in units of the reduced Planck constant

The remaining options describe the matter state of the nuclei with the given resonances.

 * `VELOCITY_DISTRIBUTION, {PARAMETERS}` (string, ...): The velocity distribution model of the nuclei in their matter state. Determines the doppler broadening of the Breit-Wigner cross section.
 * `NUCLEAR MASS` (double OR string): Mass of the nucleus in atomic mass units. It is possible to enter the numerical value or a unique identifier like '6Li'. If the identifier is entered, the value will be taken from the database of the National Institute of Standards and Technology (NIST) [[1]](#ref-nist-masses).
 * `MASS_ATTENUATION, {PARAMETERS}` (string, ...): The model for the non-resonant attenuation of photons in matter. The most convenient choice is to use the mass attenuation data from NIST [[2]](#ref-nist-attenuation) by selecting `MASS_ATTENUATION == 'nist'` and giving the file name in `mass_attenuation/` as a string parameter. The unit of the mass attenuation coefficient is `fm^2/atom` 'femtometer squared per atom'.
 * `TARGET_THICKNESS` (double): Thickness of the target in the microscopic unit `atoms/fm^2` 'atoms per femtometer squared'
 * `VELOCITY` (double): Velocity component of the target along the z-axis. This can be used to artificially doppler-shift the resonance energies.

The example `6Li_NRF` describes the single isolated resonance of 6Li at about 3563 keV. For the velocity distribution, a Maxwellian distribution with an effective temperature is chosen. The nuclear mass and non-resonant attenuation are automatically read. `6Li_RSA` shows how to stack several targets.

```
# 0) Identifier of the target (string)
6Li_NRF_scatterer
# 1) Resonance energies in eV (comma separated values)
3.562880e6
# 2) Transition widths to the ground state in eV (comma separated values)
8.16
# 3) Total transition widths in eV (comma separated values)
8.16
# 4) Angular momentum of the ground state
1
# 5) Angular momenta of the excited states (comma separated values)
0
# 6) Model for velocity distribution
maxwell_boltzmann, 412.5
# 7) Nuclear mass in atomic mass units
6Li
# 8) Mass attenuation coefficient
nist, Li.dat
# 9) Target thickness in atoms/fm^2
0.000282
# 10) Velocity of the target along the beam direction in m/s
0.
```

#### 3.3 Input file options <a name="input_file_options"></a>

TODO

### 4. License <a name="license"></a>

Copyright (C) 2018

U. Gayer (gayer.udo@gmail.com)

This code is distributed under the terms of the GNU General Public License. See the COPYING file for more information.

### 5. References <a name="references"></a>

<a name="ref-nist-masses">[1]</a> Coursey, J.S., Schwab, D.J., Tsai, J.J., and Dragoset, R.A. (2015), Atomic Weights and Isotopic Compositions (version 4.1). (https://physics.nist.gov/Comp)

<a name="ref-nist-attenuation">[2]</a> Hubbell, J.H. and Seltzer, S.M., Tables of X-Ray Mass Attenuation Coefficients
and Mass Energy-Absorption Coefficients from 1 keV to 20 MeV for Elements Z = 1 to 92 and 48 Additional Substances of Dosimetric Interest. (https://www.nist.gov/pml/x-ray-mass-attenuation-coefficients)
