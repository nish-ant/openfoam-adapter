#include "TKEGradient.H"

using namespace Foam;

preciceAdapter::FF::TKEGradient::TKEGradient(
    const Foam::fvMesh& mesh,
    const std::string nameK)
: k_(
    const_cast<volScalarField*>(
        &mesh.lookupObject<volScalarField>(nameK)))
{
    dataType_ = scalar;
}

void preciceAdapter::FF::TKEGradient::write(double* buffer, bool meshConnectivity, const unsigned int dim)
{
    int bufferIndex = 0;

    // For every boundary patch of the interface
    for (uint j = 0; j < patchIDs_.size(); j++)
    {
        int patchID = patchIDs_.at(j);

        // Get the TKE gradient boundary patch
        const scalarField gradientPatch((k_->boundaryFieldRef()[patchID]).snGrad());

        // For every cell of the patch
        forAll(gradientPatch, i)
        {
            // Copy the TKE gradient into the buffer
            buffer[bufferIndex++] = -gradientPatch[i];
        }
    }
}

void preciceAdapter::FF::TKEGradient::read(double* buffer, const unsigned int dim)
{
    int bufferIndex = 0;

    // For every boundary patch of the interface
    for (uint j = 0; j < patchIDs_.size(); j++)
    {
        int patchID = patchIDs_.at(j);

        // Get the TKE gradient boundary patch
        scalarField& gradientPatch = refCast<fixedGradientFvPatchScalarField>(
                k_->boundaryFieldRef()[patchID]).gradient();

        // For every cell of the patch
        forAll(gradientPatch, i)
        {
            // Set the TKE gradient as the buffer value
            gradientPatch[i] = buffer[bufferIndex++];
        }
    }
}

bool preciceAdapter::FF::TKEGradient::isLocationTypeSupported(const bool meshConnectivity) const
{
    return (this->locationType_ == LocationType::faceCenters);
}

std::string preciceAdapter::FF::TKEGradient::getDataName() const
{
    return "TKEGradient";
}
