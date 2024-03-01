#ifndef FF_TURBULENTSPECIFICDISSIPATIONRATE_H
#define FF_TURBULENTSPECIFICDISSIPATIONRATE_H

#include "CouplingDataUser.H"

#include "fvCFD.H"

namespace preciceAdapter
{
namespace FF
{

//- Class that writes and reads Turbulent Specific Dissipation Rate
class TurbulentSpecificDissipationRate : public CouplingDataUser
{

private:
    //- Turbulent Specific Dissipation Rate field
    Foam::volScalarField* omega_;

public:
    //- Constructor
    TurbulentSpecificDissipationRate(
        const Foam::fvMesh& mesh,
        const std::string nameOmega);

    //- Write the Turbulent Specific Dissipation Rate values into the buffer
    void write(double* buffer, bool meshConnectivity, const unsigned int dim);

    //- Read the Turbulent Specific Dissipation Rate values from the buffer
    void read(double* buffer, const unsigned int dim);

    bool isLocationTypeSupported(const bool meshConnectivity) const final;

    //- Get the name of the current data field
    std::string getDataName() const final;
};

}
}

#endif
