#ifndef FF_PRESSURERGH_GRADIENT_H
#define FF_PRESSURERGH_GRADIENT_H

#include "CouplingDataUser.H"

#include "fvCFD.H"

namespace preciceAdapter
{
namespace FF
{

//- Class that writes and reads pressureRgh gradient
class PressureRghGradient : public CouplingDataUser
{

private:
    //- PressureRgh field
    Foam::volScalarField* pRgh_;

public:
    //- Constructor
    PressureRghGradient(
        const Foam::fvMesh& mesh,
        const std::string namePRgh);

    //- Write the pressureRgh gradient values into the buffer
    void write(double* buffer, bool meshConnectivity, const unsigned int dim);

    //- Read the pressureRgh gradient values from the buffer
    void read(double* buffer, const unsigned int dim);

    bool isLocationTypeSupported(const bool meshConnectivity) const final;

    //- Get the name of the current data field
    std::string getDataName() const final;
};

}
}

#endif
