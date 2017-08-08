#ifndef CONFIG_H
#define CONFIG_H 1

#include <string>

using std::string;

//************************************
// Definitions for input/output files
//************************************

// Symbol that, when it appears at position 0 in a line, indicates a comment in all input and output files
#define COMMENT "#"

// Define directories
#define V_DIR "velocity_distribution/"
#define MU_DIR "mass_attenuation/"
#define BEAM_DIR "beam/"

// Separator between columns in a NIST x-ray attenuation file
const string NIST_SEPARATOR = "  ";

//************************************
// Definitions of input parameters
//************************************

// Number of parameters required for a single calculation
#define NPARAMETERS 10

// Length of the header of the input file, where EMIN and EMAX are set
#define INPUT_HEADER 1
// Define shortcuts to access the input parameters
#define NAME 1
#define E0 2
#define GAMMA0 3
#define GAMMA 4
#define J0 5
#define JJ 6
#define V 7
#define M 8
#define MU 9
#define Z 10

//************************************
// Definition of parameters for the calculation
//************************************

// Number of bins for the calculation
#define NBINS 1000
#define NBINS_Z 1000

// Limit for the applicability of the approximation for the doppler-shifted cross section (Gamma << Delta). The code will output a warning if the ratio Gamma/Delta is larger than APPROXIMATION_LIMIT.
#define APPROXIMATION_LIMIT 0.1 // Arbitrary choice

//************************************
// Definitions of physical constants
//************************************

#define HBARC 197.3269788e6 // in eVfm
#define HBARC2 3.89379366e16 // in eV^2fm^2
#define C 299792458. // in m/s
#define kB 8.6173303e-5 // in eV/K
#define NA 6.022140857e23 // in particles/mole
#define AtomicMassUnit 931.494095e6 // in eV/c^2
#define AtomicMassUnitG 1.660539040e-27 // in kg
#define PI 3.141592653589793

//************************************
// Parameters for plots
//************************************

#define CS_PLOT_LEGEND_X1 0.1
#define CS_PLOT_LEGEND_X2 0.4
#define CS_PLOT_LEGEND_Y1 0.8
#define CS_PLOT_LEGEND_Y2 0.9

#define V_PLOT_LEGEND_X1 0.6
#define V_PLOT_LEGEND_X2 0.9
#define V_PLOT_LEGEND_Y1 0.8
#define V_PLOT_LEGEND_Y2 0.9

#define MU_PLOT_LEGEND_X1 0.6
#define MU_PLOT_LEGEND_X2 0.9
#define MU_PLOT_LEGEND_Y1 0.8
#define MU_PLOT_LEGEND_Y2 0.9

// Maximum number of bins on the energy- und z-axis for plots of phi and alpha
#define NBINS_PHI_MAX 100
#define NBINS_Z_PHI_MAX 100

#define NBINS_ALPHA_MAX 100
#define NBINS_Z_ALPHA_MAX 100

#endif
