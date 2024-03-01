#ifndef FF_KINEMATICTHERMALCONDUCTIVITY_GRADIENT_H
#define FF_KINEMATICTHERMALCONDUCTIVITY_GRADIENT_H

#include "CouplingDataUser.H"

#include "fvCFD.H"

namespace preciceAdapter
{
namespace FF
{

//- Class that writes and reads kinematic thermal conductivity gradient
class KinematicThermalConductivityGradient : public CouplingDataUser
{

private:
    //- Kinematic thermal conductivity field
    Foam::volScalarField* kappat_;

public:
    //- Constructor
    KinematicThermalConductivityGradient(
        const Foam::fvMesh& mesh,
        const std::string nameKappat);

    //- Write the kinematic thermal conductivity gradient values into the buffer
    void write(double* buffer, bool meshConnectivity, const unsigned int dim);

    //- Read the kinematic thermal conductivity gradient values from the buffer
    void read(double* buffer, const unsigned int dim);

    bool isLocationTypeSupported(const bool meshConnectivity) const final;

    //- Get the name of the current data field
    std::string getDataName() const final;
};

}
}

#endif
