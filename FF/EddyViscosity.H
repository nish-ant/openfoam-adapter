#ifndef FF_EDDYVISCOSITY_H
#define FF_EDDYVISCOSITY_H

#include "CouplingDataUser.H"

#include "fvCFD.H"

namespace preciceAdapter
{
namespace FF
{

//- Class that writes and reads eddy viscosity
class EddyViscosity : public CouplingDataUser
{

private:
    //- Eddy viscosity field
    Foam::volScalarField* nut_;

public:
    //- Constructor
    EddyViscosity(
        const Foam::fvMesh& mesh,
        const std::string nameNut);

    //- Write the eddy viscosity values into the buffer
    void write(double* buffer, bool meshConnectivity, const unsigned int dim);

    //- Read the eddy viscosity values from the buffer
    void read(double* buffer, const unsigned int dim);

    bool isLocationTypeSupported(const bool meshConnectivity) const final;

    //- Get the name of the current data field
    std::string getDataName() const final;
};

}
}

#endif
