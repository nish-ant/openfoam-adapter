#include "TurbulentSpecificDissipationRate.H"

using namespace Foam;

preciceAdapter::FF::TurbulentSpecificDissipationRate::TurbulentSpecificDissipationRate(
    const Foam::fvMesh& mesh,
    const std::string nameOmega)
: omega_(
    const_cast<volScalarField*>(
        &mesh.lookupObject<volScalarField>(nameOmega)))
{
    dataType_ = scalar;
}

void preciceAdapter::FF::TurbulentSpecificDissipationRate::write(double* buffer, bool meshConnectivity, const unsigned int dim)
{
    int bufferIndex = 0;

    // For every boundary patch of the interface
    for (uint j = 0; j < patchIDs_.size(); j++)
    {
        int patchID = patchIDs_.at(j);

        // For every cell of the patch
        forAll(omega_->boundaryFieldRef()[patchID], i)
        {
            // Copy the Turbulent Specific Dissipation Rate into the buffer
            buffer[bufferIndex++] = omega_->boundaryFieldRef()[patchID][i];
        }
    }
}

void preciceAdapter::FF::TurbulentSpecificDissipationRate::read(double* buffer, const unsigned int dim)
{
    int bufferIndex = 0;

    // For every boundary patch of the interface
    for (uint j = 0; j < patchIDs_.size(); j++)
    {
        int patchID = patchIDs_.at(j);

        // For every cell of the patch
        forAll(omega_->boundaryFieldRef()[patchID], i)
        {
            // Set the Turbulent Specific Dissipation Rate as the buffer value
            omega_->boundaryFieldRef()[patchID][i] = buffer[bufferIndex++];
        }
    }
}

bool preciceAdapter::FF::TurbulentSpecificDissipationRate::isLocationTypeSupported(const bool meshConnectivity) const
{
    return (this->locationType_ == LocationType::faceCenters);
}

std::string preciceAdapter::FF::TurbulentSpecificDissipationRate::getDataName() const
{
    return "TurbulentSpecificDissipationRate";
}
