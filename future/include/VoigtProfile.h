#pragma once

/// \brief Class for the Voigt profile.
///
/// The Voigt profile \cite VoigtProfile2019 is a convolution of a Breit-Wigner distribution
/// and a normal distribution. It describes the lineshape of a Doppler-broadened nuclear resonance
/// in the case that the velocities of the nuclei are normal distributed.
class VoigtProfile{

public:
    /// \brief Default constructor.
    VoigtProfile() = default;
    /// \brief Default destructor.
    ~VoigtProfile() = default;

    /// \brief Function to evaluate the probability distribution function (PDF).
    /// \param x Point at which the PDF should be evaluated.
    /// \param sigma Standard deviation of the normal distribution.
    /// \param gamma Width of the Breit-Wigner distribution.
    ///
    /// At the moment, the PDF of the Voigt profile is taken from the 
    /// ROOT data analysis framework \cite Brun1996 \cite Brun1997.
    double pdf( const double x, const double sigma, const double gamma ) const;

    /// \brief Function to evaluate the cumulative distribution function (CDF).
    ///
    /// \todo Find an implementation or create own implementation. At the moment, the grids for the evaluation of a Voigt-type distribution are generated either for a Breit-Wigner or a normal distribution, depending on the relative magnitude of `sigma` and `gamma`. A more general method, based on the actual shape of the distribution, would be desirable.
    double cdf( const double x, const double sigma, const double gamma ) const;
};