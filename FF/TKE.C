#include "TKE.H"

using namespace Foam;

preciceAdapter::FF::TKE::TKE(
    const Foam::fvMesh& mesh,
    const std::string nameK)
: k_(
    const_cast<volScalarField*>(
        &mesh.lookupObject<volScalarField>(nameK)))
{
    dataType_ = scalar;
}

void preciceAdapter::FF::TKE::write(double* buffer, bool meshConnectivity, const unsigned int dim)
{
    int bufferIndex = 0;

    // For every boundary patch of the interface
    for (uint j = 0; j < patchIDs_.size(); j++)
    {
        int patchID = patchIDs_.at(j);

        // For every cell of the patch
        forAll(k_->boundaryFieldRef()[patchID], i)
        {
            // Copy the TKE into the buffer
            buffer[bufferIndex++] = k_->boundaryFieldRef()[patchID][i];
        }
    }
}

void preciceAdapter::FF::TKE::read(double* buffer, const unsigned int dim)
{
    int bufferIndex = 0;

    // For every boundary patch of the interface
    for (uint j = 0; j < patchIDs_.size(); j++)
    {
        int patchID = patchIDs_.at(j);

        // For every cell of the patch
        forAll(k_->boundaryFieldRef()[patchID], i)
        {
            // Set the TKE as the buffer value
            k_->boundaryFieldRef()[patchID][i] = buffer[bufferIndex++];
        }
    }
}

bool preciceAdapter::FF::TKE::isLocationTypeSupported(const bool meshConnectivity) const
{
    return (this->locationType_ == LocationType::faceCenters);
}

std::string preciceAdapter::FF::TKE::getDataName() const
{
    return "TKE";
}
