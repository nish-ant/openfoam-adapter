#ifndef FF_EDDYVISCOSITY_GRADIENT_H
#define FF_EDDYVISCOSITY_GRADIENT_H

#include "CouplingDataUser.H"

#include "fvCFD.H"

namespace preciceAdapter
{
namespace FF
{

//- Class that writes and reads eddy viscosity gradient
class EddyViscosityGradient : public CouplingDataUser
{

private:
    //- Eddy viscosity field
    Foam::volScalarField* nut_;

public:
    //- Constructor
    EddyViscosityGradient(
        const Foam::fvMesh& mesh,
        const std::string nameNut);

    //- Write the eddy viscosity gradient values into the buffer
    void write(double* buffer, bool meshConnectivity, const unsigned int dim);

    //- Read the eddy viscosity gradient values from the buffer
    void read(double* buffer, const unsigned int dim);

    bool isLocationTypeSupported(const bool meshConnectivity) const final;

    //- Get the name of the current data field
    std::string getDataName() const final;
};

}
}

#endif
