#include "PressureRgh.H"

using namespace Foam;

preciceAdapter::FF::PressureRgh::PressureRgh(
    const Foam::fvMesh& mesh,
    const std::string namePRgh)
: pRgh_(
    const_cast<volScalarField*>(
        &mesh.lookupObject<volScalarField>(namePRgh)))
{
    dataType_ = scalar;
}

void preciceAdapter::FF::PressureRgh::write(double* buffer, bool meshConnectivity, const unsigned int dim)
{
    int bufferIndex = 0;

    // For every boundary patch of the interface
    for (uint j = 0; j < patchIDs_.size(); j++)
    {
        int patchID = patchIDs_.at(j);

        // For every cell of the patch
        forAll(pRgh_->boundaryFieldRef()[patchID], i)
        {
            // Copy the pressureRgh into the buffer
            buffer[bufferIndex++] = pRgh_->boundaryFieldRef()[patchID][i];
        }
    }
}

void preciceAdapter::FF::PressureRgh::read(double* buffer, const unsigned int dim)
{
    int bufferIndex = 0;

    // For every boundary patch of the interface
    for (uint j = 0; j < patchIDs_.size(); j++)
    {
        int patchID = patchIDs_.at(j);

        // For every cell of the patch
        forAll(pRgh_->boundaryFieldRef()[patchID], i)
        {
            // Set the pressureRgh as the buffer value
            pRgh_->boundaryFieldRef()[patchID][i] = buffer[bufferIndex++];
        }
    }
}

bool preciceAdapter::FF::PressureRgh::isLocationTypeSupported(const bool meshConnectivity) const
{
    return (this->locationType_ == LocationType::faceCenters);
}

std::string preciceAdapter::FF::PressureRgh::getDataName() const
{
    return "PressureRgh";
}
