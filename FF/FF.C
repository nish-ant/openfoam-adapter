#include "FF.H"

#include "Utilities.H"

using namespace Foam;

preciceAdapter::FF::FluidFluid::FluidFluid(
    const Foam::fvMesh& mesh)
: mesh_(mesh)
{
}

bool preciceAdapter::FF::FluidFluid::configure(const IOdictionary& adapterConfig)
{
    DEBUG(adapterInfo("Configuring the FF module..."));

    // Read the FF-specific options from the adapter's configuration file
    if (!readConfig(adapterConfig))
    {
        return false;
    }

    // NOTE: If you want to add a new solver type, which you can manually
    // specify in the configuration, add it here. See also the methods
    // addWriters() and addReaders().
    // Check the solver type and determine it if needed
    if (
        solverType_.compare("compressible") == 0 || solverType_.compare("incompressible") == 0)
    {
        DEBUG(adapterInfo("Known solver type: " + solverType_));
    }
    else if (solverType_.compare("none") == 0)
    {
        DEBUG(adapterInfo("Determining the solver type..."));
        solverType_ = determineSolverType();
    }
    else
    {
        DEBUG(adapterInfo("Determining the solver type for the FF module... (override by setting solverType to one of {compressible, incompressible})"));
        solverType_ = determineSolverType();
    }

    return true;
}

bool preciceAdapter::FF::FluidFluid::readConfig(const IOdictionary& adapterConfig)
{
    const dictionary& FFdict = adapterConfig.subOrEmptyDict("FF");

    // Read the solver type (if not specified, it is determined automatically)
    solverType_ = FFdict.lookupOrDefault<word>("solverType", "");
    DEBUG(adapterInfo("    user-defined solver type : " + solverType_));

    // Read the name of the velocity field (if different)
    nameU_ = FFdict.lookupOrDefault<word>("nameU", "U");
    DEBUG(adapterInfo("    velocity field name : " + nameU_));

    // Read the name of the pressure field (if different)
    nameP_ = FFdict.lookupOrDefault<word>("nameP", "p");
    DEBUG(adapterInfo("    pressure field name : " + nameP_));

    // Read the name of the pressureRgh field (if different)
    namePRgh_ = FFdict.lookupOrDefault<word>("namePRgh", "p_rgh");
    DEBUG(adapterInfo("    pressureRgh field name : " + namePRgh_));

    // Read the name of the temperature field (if different)
    nameT_ = FFdict.lookupOrDefault<word>("nameT", "T");
    DEBUG(adapterInfo("    temperature field name : " + nameT_));

    // Read the name of the eddy viscosity field (if different)
    nameNut_ = FFdict.lookupOrDefault<word>("nameNut", "nut");
    DEBUG(adapterInfo("    eddy viscosity field name : " + nameNut_));

    // Read the name of the kinematic thermal conductivity field (if different)
    nameKappat_ = FFdict.lookupOrDefault<word>("nameKappat", "kappat");
    DEBUG(adapterInfo("    kinematic thermal conductivity field name : " + nameKappat_));

    // Read the name of the turbulent kinetic energy field (if different)
    nameK_ = FFdict.lookupOrDefault<word>("nameK", "k");
    DEBUG(adapterInfo("    turbulent kinetic energy field name : " + nameK_));

    // Read the name of the turbulent dissipation rate field (if different)
    nameEpsilon_ = FFdict.lookupOrDefault<word>("nameEpsilon", "epsilon");
    DEBUG(adapterInfo("    turbulent dissipation rate field name : " + nameEpsilon_));

    // Read the name of the specific turbulent dissipation rate field (if different)
    nameOmega_ = FFdict.lookupOrDefault<word>("nameOmega", "omega");
    DEBUG(adapterInfo("    turbulent specific dissipation rate field name : " + nameOmega_));

    return true;
}

std::string preciceAdapter::FF::FluidFluid::determineSolverType()
{
    // NOTE: When coupling a different variable, you may want to
    // add more cases here. Or you may provide the solverType in the config.

    std::string solverType = "unknown";

    // Dimensions: pressure or pressureRgh
    dimensionSet dimensionsCompressible(1, -1, -2, 0, 0, 0, 0);
    dimensionSet dimensionsIncompressible(0, 2, -2, 0, 0, 0, 0);

    if (mesh_.foundObject<volScalarField>("p"))
    {
        const volScalarField& p_ = mesh_.lookupObject<volScalarField>("p");

        if (p_.dimensions() == dimensionsCompressible)
            solverType = "compressible";
        else if (p_.dimensions() == dimensionsIncompressible)
            solverType = "incompressible";
    }

    else if (mesh_.foundObject<volScalarField>("p_rgh"))
    {
        const volScalarField& pRgh_ = mesh_.lookupObject<volScalarField>("p_rgh");

        if (pRgh_.dimensions() == dimensionsCompressible)
            solverType = "compressible";
        else if (pRgh_.dimensions() == dimensionsIncompressible)
            solverType = "incompressible";
    }

    // TODO: Add special case for multiphase solvers.
    // Currently, interFoam is misclassified as "compressible".

    if (solverType == "unknown")
        adapterInfo("Failed to determine the solver type. "
                    "Please specify your solver type in the FF section of the "
                    "preciceDict. Known solver types for FF are: "
                    "incompressible and "
                    "compressible",
                    "error");

    DEBUG(adapterInfo("Automatically determined solver type : " + solverType));

    return solverType;
}

bool preciceAdapter::FF::FluidFluid::addWriters(std::string dataName, Interface* interface)
{
    bool found = true; // Set to false later, if needed.

    if (dataName.find("VelocityGradient") == 0)
    {
        interface->addCouplingDataWriter(
            dataName,
            new VelocityGradient(mesh_, nameU_));
        DEBUG(adapterInfo("Added writer: Velocity Gradient."));
    }
    else if (dataName.find("Velocity") == 0)
    {
        interface->addCouplingDataWriter(
            dataName,
            new Velocity(mesh_, nameU_));
        DEBUG(adapterInfo("Added writer: Velocity."));
    }

    else if (dataName.find("PressureGradient") == 0)
    {
        interface->addCouplingDataWriter(
            dataName,
            new PressureGradient(mesh_, nameP_));
        DEBUG(adapterInfo("Added writer: Pressure Gradient."));
    }
    else if (dataName.find("Pressure") == 0)
    {
        interface->addCouplingDataWriter(
            dataName,
            new Pressure(mesh_, nameP_));
        DEBUG(adapterInfo("Added writer: Pressure."));
    }

    else if (dataName.find("PressureRghGradient") == 0)
    {
        interface->addCouplingDataWriter(
            dataName,
            new PressureRghGradient(mesh_, namePRgh_));
        DEBUG(adapterInfo("Added writer: PressureRgh Gradient."));
    }
    else if (dataName.find("PressureRgh") == 0)
    {
        interface->addCouplingDataWriter(
            dataName,
            new PressureRgh(mesh_, namePRgh_));
        DEBUG(adapterInfo("Added writer: PressureRgh."));
    }

    else if (dataName.find("TemperatureGradient") == 0)
    {
        interface->addCouplingDataWriter(
            dataName,
            new TemperatureGradient(mesh_, nameT_));
        DEBUG(adapterInfo("Added writer: Temperature Gradient."));
    }
    else if (dataName.find("Temperature") == 0)
    {
        interface->addCouplingDataWriter(
            dataName,
            new Temperature(mesh_, nameT_));
        DEBUG(adapterInfo("Added writer: Temperature."));
    }

    else if (dataName.find("EddyViscosityGradient") == 0)
    {
        interface->addCouplingDataWriter(
            dataName,
            new EddyViscosityGradient(mesh_, nameNut_));
        DEBUG(adapterInfo("Added writer: Eddy Viscosity Gradient."));
    }
    else if (dataName.find("EddyViscosity") == 0)
    {
        interface->addCouplingDataWriter(
            dataName,
            new EddyViscosity(mesh_, nameNut_));
        DEBUG(adapterInfo("Added writer: Eddy Viscosity."));
    }

    else if (dataName.find("KinematicThermalConductivityGradient") == 0)
    {
        interface->addCouplingDataWriter(
            dataName,
            new KinematicThermalConductivityGradient(mesh_, nameKappat_));
        DEBUG(adapterInfo("Added writer: Kinematic Thermal Conductivity Gradient."));
    }
    else if (dataName.find("KinematicThermalConductivity") == 0)
    {
        interface->addCouplingDataWriter(
            dataName,
            new KinematicThermalConductivity(mesh_, nameKappat_));
        DEBUG(adapterInfo("Added writer: Kinematic Thermal Conductivity."));
    }

    else if (dataName.find("TKEGradient") == 0)
    {
        interface->addCouplingDataWriter(
            dataName,
            new TKEGradient(mesh_, nameK_));
        DEBUG(adapterInfo("Added writer: TKE Gradient."));
    }
    else if (dataName.find("TKE") == 0)
    {
        interface->addCouplingDataWriter(
            dataName,
            new TKE(mesh_, nameK_));
        DEBUG(adapterInfo("Added writer: TKE."));
    }

    else if (dataName.find("TurbulentDissipationRateGradient") == 0)
    {
        interface->addCouplingDataWriter(
            dataName,
            new TurbulentDissipationRateGradient(mesh_, nameEpsilon_));
        DEBUG(adapterInfo("Added writer: Turbulent Dissipation Rate Gradient."));
    }
    else if (dataName.find("TurbulentDissipationRate") == 0)
    {
        interface->addCouplingDataWriter(
            dataName,
            new TurbulentDissipationRate(mesh_, nameEpsilon_));
        DEBUG(adapterInfo("Added writer: Turbulent Dissipation Rate."));
    }

    else if (dataName.find("TurbulentSpecificDissipationRateGradient") == 0)
    {
        interface->addCouplingDataWriter(
            dataName,
            new TurbulentSpecificDissipationRateGradient(mesh_, nameOmega_));
        DEBUG(adapterInfo("Added writer: Turbulent Specific Dissipation Rate Gradient."));
    }
    else if (dataName.find("TurbulentSpecificDissipationRate") == 0)
    {
        interface->addCouplingDataWriter(
            dataName,
            new TurbulentSpecificDissipationRate(mesh_, nameOmega_));
        DEBUG(adapterInfo("Added writer: Turbulent Specific Dissipation Rate."));
    }

    else
    {
        found = false;
    }

    // NOTE: If you want to couple another variable, you need
    // to add your new coupling data user as a coupling data
    // writer here (and as a reader below).
    // The argument of the dataName.compare() needs to match
    // the one provided in the adapter's configuration file.

    return found;
}

bool preciceAdapter::FF::FluidFluid::addReaders(std::string dataName, Interface* interface)
{
    bool found = true; // Set to false later, if needed.

    if (dataName.find("VelocityGradient") == 0)
    {
        interface->addCouplingDataReader(
            dataName,
            new VelocityGradient(mesh_, nameU_));
        DEBUG(adapterInfo("Added reader: VelocityGradient."));
    }
    else if (dataName.find("Velocity") == 0)
    {
        interface->addCouplingDataReader(
            dataName,
            new Velocity(mesh_, nameU_));
        DEBUG(adapterInfo("Added reader: Velocity."));
    }

    else if (dataName.find("PressureGradient") == 0)
    {
        interface->addCouplingDataReader(
            dataName,
            new PressureGradient(mesh_, nameP_));
        DEBUG(adapterInfo("Added reader: Pressure Gradient."));
    }
    else if (dataName.find("Pressure") == 0)
    {
        interface->addCouplingDataReader(
            dataName,
            new Pressure(mesh_, nameP_));
        DEBUG(adapterInfo("Added reader: Pressure."));
    }

    else if (dataName.find("PressureRghGradient") == 0)
    {
        interface->addCouplingDataReader(
            dataName,
            new PressureRghGradient(mesh_, namePRgh_));
        DEBUG(adapterInfo("Added reader: PressureRgh Gradient."));
    }
    else if (dataName.find("PressureRgh") == 0)
    {
        interface->addCouplingDataReader(
            dataName,
            new PressureRgh(mesh_, namePRgh_));
        DEBUG(adapterInfo("Added reader: PressureRgh."));
    }

    else if (dataName.find("TemperatureGradient") == 0)
    {
        interface->addCouplingDataReader(
            dataName,
            new TemperatureGradient(mesh_, nameT_));
        DEBUG(adapterInfo("Added reader: Temperature Gradient."));
    }
    else if (dataName.find("Temperature") == 0)
    {
        interface->addCouplingDataReader(
            dataName,
            new Temperature(mesh_, nameT_));
        DEBUG(adapterInfo("Added reader: Temperature."));
    }

    else if (dataName.find("EddyViscosityGradient") == 0)
    {
        interface->addCouplingDataReader(
            dataName,
            new EddyViscosityGradient(mesh_, nameNut_));
        DEBUG(adapterInfo("Added reader: Eddy Viscosity Gradient."));
    }
    else if (dataName.find("EddyViscosity") == 0)
    {
        interface->addCouplingDataReader(
            dataName,
            new EddyViscosity(mesh_, nameNut_));
        DEBUG(adapterInfo("Added reader: Eddy Viscosity."));
    }

    else if (dataName.find("KinematicThermalConductivityGradient") == 0)
    {
        interface->addCouplingDataReader(
            dataName,
            new KinematicThermalConductivityGradient(mesh_, nameKappat_));
        DEBUG(adapterInfo("Added reader: Kinematic Thermal Conductivity Gradient."));
    }
    else if (dataName.find("KinematicThermalConductivity") == 0)
    {
        interface->addCouplingDataReader(
            dataName,
            new KinematicThermalConductivity(mesh_, nameKappat_));
        DEBUG(adapterInfo("Added reader: Kinematic Thermal Conductivity."));
    }

    else if (dataName.find("TKEGradient") == 0)
    {
        interface->addCouplingDataReader(
            dataName,
            new TKEGradient(mesh_, nameK_));
        DEBUG(adapterInfo("Added reader: TKE Gradient."));
    }
    else if (dataName.find("TKE") == 0)
    {
        interface->addCouplingDataReader(
            dataName,
            new TKE(mesh_, nameK_));
        DEBUG(adapterInfo("Added reader: TKE."));
    }

    else if (dataName.find("TurbulentDissipationRateGradient") == 0)
    {
        interface->addCouplingDataReader(
            dataName,
            new TurbulentDissipationRateGradient(mesh_, nameEpsilon_));
        DEBUG(adapterInfo("Added reader: Turbulent Dissipation Rate Gradient."));
    }
    else if (dataName.find("TurbulentDissipationRate") == 0)
    {
        interface->addCouplingDataReader(
            dataName,
            new TurbulentDissipationRate(mesh_, nameEpsilon_));
        DEBUG(adapterInfo("Added reader: Turbulent Dissipation Rate."));
    }

    else if (dataName.find("TurbulentSpecificDissipationRateGradient") == 0)
    {
        interface->addCouplingDataReader(
            dataName,
            new TurbulentSpecificDissipationRateGradient(mesh_, nameOmega_));
        DEBUG(adapterInfo("Added reader: Turbulent Specific Dissipation Rate Gradient."));
    }
    else if (dataName.find("TurbulentSpecificDissipationRate") == 0)
    {
        interface->addCouplingDataReader(
            dataName,
            new TurbulentSpecificDissipationRate(mesh_, nameOmega_));
        DEBUG(adapterInfo("Added reader: Turbulent Specific Dissipation Rate."));
    }

    else
    {
        found = false;
    }

    // NOTE: If you want to couple another variable, you need
    // to add your new coupling data user as a coupling data
    // reader here (and as a writer above).
    // The argument of the dataName.compare() needs to match
    // the one provided in the adapter's configuration file.

    return found;
}
