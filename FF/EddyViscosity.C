#include "EddyViscosity.H"

using namespace Foam;

preciceAdapter::FF::EddyViscosity::EddyViscosity(
    const Foam::fvMesh& mesh,
    const std::string nameNut)
: nut_(
    const_cast<volScalarField*>(
        &mesh.lookupObject<volScalarField>(nameNut)))
{
    dataType_ = scalar;
}

void preciceAdapter::FF::EddyViscosity::write(double* buffer, bool meshConnectivity, const unsigned int dim)
{
    int bufferIndex = 0;

    // For every boundary patch of the interface
    for (uint j = 0; j < patchIDs_.size(); j++)
    {
        int patchID = patchIDs_.at(j);

        // For every cell of the patch
        forAll(nut_->boundaryFieldRef()[patchID], i)
        {
            // Copy the eddy viscosity into the buffer
            buffer[bufferIndex++] = nut_->boundaryFieldRef()[patchID][i];
        }
    }
}

void preciceAdapter::FF::EddyViscosity::read(double* buffer, const unsigned int dim)
{
    int bufferIndex = 0;

    // For every boundary patch of the interface
    for (uint j = 0; j < patchIDs_.size(); j++)
    {
        int patchID = patchIDs_.at(j);

        // For every cell of the patch
        forAll(nut_->boundaryFieldRef()[patchID], i)
        {
            // Set the eddy viscosity as the buffer value
            nut_->boundaryFieldRef()[patchID][i] = buffer[bufferIndex++];
        }
    }
}

bool preciceAdapter::FF::EddyViscosity::isLocationTypeSupported(const bool meshConnectivity) const
{
    return (this->locationType_ == LocationType::faceCenters);
}

std::string preciceAdapter::FF::EddyViscosity::getDataName() const
{
    return "EddyViscosity";
}
