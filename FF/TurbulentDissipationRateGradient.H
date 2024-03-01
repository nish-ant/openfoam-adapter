#ifndef FF_TURBULENTDISSIPATIONRATE_GRADIENT_H
#define FF_TURBULENTDISSIPATIONRATE_GRADIENT_H

#include "CouplingDataUser.H"

#include "fvCFD.H"

namespace preciceAdapter
{
namespace FF
{

//- Class that writes and reads Turbulent Dissipation Rate gradient
class TurbulentDissipationRateGradient : public CouplingDataUser
{

private:
    //- Turbulent Dissipation Rate field
    Foam::volScalarField* epsilon_;

public:
    //- Constructor
    TurbulentDissipationRateGradient(
        const Foam::fvMesh& mesh,
        const std::string nameEpsilon);

    //- Write the Turbulent Dissipation Rate gradient values into the buffer
    void write(double* buffer, bool meshConnectivity, const unsigned int dim);

    //- Read the Turbulent Dissipation Rate gradient values from the buffer
    void read(double* buffer, const unsigned int dim);

    bool isLocationTypeSupported(const bool meshConnectivity) const final;

    //- Get the name of the current data field
    std::string getDataName() const final;
};

}
}

#endif
