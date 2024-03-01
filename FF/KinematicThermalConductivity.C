#include "KinematicThermalConductivity.H"

using namespace Foam;

preciceAdapter::FF::KinematicThermalConductivity::KinematicThermalConductivity(
    const Foam::fvMesh& mesh,
    const std::string nameKappat)
: kappat_(
    const_cast<volScalarField*>(
        &mesh.lookupObject<volScalarField>(nameKappat)))
{
    dataType_ = scalar;
}

void preciceAdapter::FF::KinematicThermalConductivity::write(double* buffer, bool meshConnectivity, const unsigned int dim)
{
    int bufferIndex = 0;

    // For every boundary patch of the interface
    for (uint j = 0; j < patchIDs_.size(); j++)
    {
        int patchID = patchIDs_.at(j);

        // For every cell of the patch
        forAll(kappat_->boundaryFieldRef()[patchID], i)
        {
            // Copy the kinematic thermal conductivity into the buffer
            buffer[bufferIndex++] = kappat_->boundaryFieldRef()[patchID][i];
        }
    }
}

void preciceAdapter::FF::KinematicThermalConductivity::read(double* buffer, const unsigned int dim)
{
    int bufferIndex = 0;

    // For every boundary patch of the interface
    for (uint j = 0; j < patchIDs_.size(); j++)
    {
        int patchID = patchIDs_.at(j);

        // For every cell of the patch
        forAll(kappat_->boundaryFieldRef()[patchID], i)
        {
            // Set the kinematic thermal conductivity as the buffer value
            kappat_->boundaryFieldRef()[patchID][i] = buffer[bufferIndex++];
        }
    }
}

bool preciceAdapter::FF::KinematicThermalConductivity::isLocationTypeSupported(const bool meshConnectivity) const
{
    return (this->locationType_ == LocationType::faceCenters);
}

std::string preciceAdapter::FF::KinematicThermalConductivity::getDataName() const
{
    return "KinematicThermalConductivity";
}
