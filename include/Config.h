/*    
    This file is part of SeAN.

    SeAN is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    SeAN is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with SeAN.  If not, see <http://www.gnu.org/licenses/>.
*/


#ifndef CONFIG_H
#define CONFIG_H 1

#include <string>

using std::string;

// Definitions of physical constants
const double HBARC = 197.3269788e6; // in eVfm
const double HBARC2 = 3.89379366e16; // in eV^2fm^2
const double SPEEDOFLIGHT = 299792458.; // in m/s
const double kB = 8.6173303e-5; // in eV/K
const double AtomicMassUnit = 931.494095e6; // in eV/c^2
const double AtomicMassUnitG = 1.660539040e-24; // in g
const double PI = 3.141592653589793;
const double PI2 = 9.869604401089358; 

// Settings for the SeAN input file
const string COMMENT = "#";
const char DELIMITER = ',';
const char LOOP_DELIMITER = ';';
const char LOOP_START = '[';
const char LOOP_STOP = ']';
const char SIM_LOOP_START = '{';
const char SIM_LOOP_STOP = '}';
const unsigned int N_TARGET_SETTINGS = 11;

// Settings for the directories

const string VELOCITY_DISTRIBUTION_DIR = "velocity_distribution/";
const string PHONON_DIR = "phdos/";
const string CROSS_SECTION_DIR = "cross_section/";
const string MU_DIR = "mass_attenuation/";
const string BEAM_DIR = "beam/";
const string MASS_DIR = "atomic_mass/";
const string AME_FILE_NAME = "mass16.txt";

// Settings for phonon calculations
const unsigned int MU_MAX_INTEGRAL = 3; // Set integration range for mu in terms of the 'decay constant' 2/Gamma 

// Settings for SeAN output files
const string HORIZONTAL_LINE = "##############################################################";
const string TXT_SUFFIX = ".txt";
const string CAL_SUFFIX = ".cal";
const string CAL_FILE_NAME = "calibration";

#endif
