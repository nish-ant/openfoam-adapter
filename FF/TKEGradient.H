#ifndef FF_TKE_GRADIENT_H
#define FF_TKE_GRADIENT_H

#include "CouplingDataUser.H"

#include "fvCFD.H"

namespace preciceAdapter
{
namespace FF
{

//- Class that writes and reads TKE gradient
class TKEGradient : public CouplingDataUser
{

private:
    //- TKE field
    Foam::volScalarField* k_;

public:
    //- Constructor
    TKEGradient(
        const Foam::fvMesh& mesh,
        const std::string nameK);

    //- Write the TKE gradient values into the buffer
    void write(double* buffer, bool meshConnectivity, const unsigned int dim);

    //- Read the TKE gradient values from the buffer
    void read(double* buffer, const unsigned int dim);

    bool isLocationTypeSupported(const bool meshConnectivity) const final;

    //- Get the name of the current data field
    std::string getDataName() const final;
};

}
}

#endif
