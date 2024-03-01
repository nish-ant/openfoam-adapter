#ifndef FF_H
#define FF_H

#include "Interface.H"

#include "FF/Velocity.H"
#include "FF/Pressure.H"
#include "FF/VelocityGradient.H"
#include "FF/PressureGradient.H"

#include "FF/PressureRgh.H"
#include "FF/Temperature.H"
#include "FF/EddyViscosity.H"
#include "FF/KinematicThermalConductivity.H"
#include "FF/TKE.H"
#include "FF/TurbulentDissipationRate.H"
#include "FF/TurbulentSpecificDissipationRate.H"

#include "FF/PressureRghGradient.H"
#include "FF/TemperatureGradient.H"
#include "FF/EddyViscosityGradient.H"
#include "FF/KinematicThermalConductivityGradient.H"
#include "FF/TKEGradient.H"
#include "FF/TurbulentDissipationRateGradient.H"
#include "FF/TurbulentSpecificDissipationRateGradient.H"

#include "fvCFD.H"

namespace preciceAdapter
{
namespace FF
{

class FluidFluid
{

protected:
    //- OpenFOAM fvMesh object
    const Foam::fvMesh& mesh_;

    // TODO: Create a struct for all the parameter names

    //- Solver type
    std::string solverType_ = "none";

    //- Name of the velocity field
    std::string nameU_ = "U";

    //- Name of the pressure field
    std::string nameP_ = "p";

    //- Name of the pressureRgh field
    std::string namePRgh_ = "p_rgh";

    //- Name of the temperature field
    std::string nameT_ = "T";

    //- Name of the eddy viscosity field
    std::string nameNut_ = "nut";

    //- Name of the kinematic thermal conductivity field
    std::string nameKappat_ = "kappat";

    //- Name of the TKE field
    std::string nameK_ = "k";

    //- Name of the turbulent dissipation rate field
    std::string nameEpsilon_ = "epsilon";

    //- Name of the turbulent specific dissipation rate field
    std::string nameOmega_ = "omega";

    //- Determine the solver type
    std::string determineSolverType();

    //- Read the FF-related options from the adapter's configuration file
    bool readConfig(const IOdictionary& adapterConfig);

public:
    //- Constructor
    FluidFluid(const Foam::fvMesh& mesh);

    //- Configure
    bool configure(const IOdictionary& adapterConfig);

    //- Add coupling data writers
    bool addWriters(std::string dataName, Interface* interface);

    //- Add coupling data readers
    bool addReaders(std::string dataName, Interface* interface);
};

}
}

#endif
