# SeAN
## Self Absorption Numerical

## Table of Contents

1. [Introduction](#introduction)
2. [Installation](#installation)
	2.1 [Prerequisites](#prerequisites)
	2.2 [Compilation](#compilation)
	2.3 [Build options](build_options)
3. [Usage](#usage)
4. [License](#license)
5. [References](#references)

### 1. Introduction <a name="introduction"></a>

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

By default, the program is built in debug mode. A small increase in computation speed can be gained by replacing the `-g` compiler flag with an optimization flag `-O1`, `-O2` or `-O3`. No extensive testing has been done by the author as to which of the three is the best tradeoff between compilation time and speedup, at the moment.

##### Change C++ compiler

If you don't have the default compiler `g++` or want to use another one on purpose, change the value of the `CC` variable in the Makefile.

### 3. Usage <a name="usage"></a>

### 4. License <a name="license"></a>

### 5. References <a name="references"></a>