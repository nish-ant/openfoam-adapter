#include "TurbulentDissipationRate.H"

using namespace Foam;

preciceAdapter::FF::TurbulentDissipationRate::TurbulentDissipationRate(
    const Foam::fvMesh& mesh,
    const std::string nameEpsilon)
: epsilon_(
    const_cast<volScalarField*>(
        &mesh.lookupObject<volScalarField>(nameEpsilon)))
{
    dataType_ = scalar;
}

void preciceAdapter::FF::TurbulentDissipationRate::write(double* buffer, bool meshConnectivity, const unsigned int dim)
{
    int bufferIndex = 0;

    // For every boundary patch of the interface
    for (uint j = 0; j < patchIDs_.size(); j++)
    {
        int patchID = patchIDs_.at(j);

        // For every cell of the patch
        forAll(epsilon_->boundaryFieldRef()[patchID], i)
        {
            // Copy the Turbulent Dissipation Rate into the buffer
            buffer[bufferIndex++] = epsilon_->boundaryFieldRef()[patchID][i];
        }
    }
}

void preciceAdapter::FF::TurbulentDissipationRate::read(double* buffer, const unsigned int dim)
{
    int bufferIndex = 0;

    // For every boundary patch of the interface
    for (uint j = 0; j < patchIDs_.size(); j++)
    {
        int patchID = patchIDs_.at(j);

        // For every cell of the patch
        forAll(epsilon_->boundaryFieldRef()[patchID], i)
        {
            // Set the Turbulent Dissipation Rate as the buffer value
            epsilon_->boundaryFieldRef()[patchID][i] = buffer[bufferIndex++];
        }
    }
}

bool preciceAdapter::FF::TurbulentDissipationRate::isLocationTypeSupported(const bool meshConnectivity) const
{
    return (this->locationType_ == LocationType::faceCenters);
}

std::string preciceAdapter::FF::TurbulentDissipationRate::getDataName() const
{
    return "TurbulentDissipationRate";
}
