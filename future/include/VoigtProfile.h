#pragma once

class VoigtProfile{

public:
    VoigtProfile() = default;
    ~VoigtProfile() = default;

    double pdf( double x, double sigma, double gamma ) const;
    double cdf( double x ) const;
};