#include "TMath.h"

#include "VoigtProfile.h"

double VoigtProfile::pdf( double x, double sigma, double gamma ) const {
    return TMath::Voigt(x, sigma, gamma);
}