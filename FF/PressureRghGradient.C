#include "PressureRghGradient.H"

using namespace Foam;

preciceAdapter::FF::PressureRghGradient::PressureRghGradient(
    const Foam::fvMesh& mesh,
    const std::string namePRgh)
: pRgh_(
    const_cast<volScalarField*>(
        &mesh.lookupObject<volScalarField>(namePRgh)))
{
    dataType_ = scalar;
}

void preciceAdapter::FF::PressureRghGradient::write(double* buffer, bool meshConnectivity, const unsigned int dim)
{
    int bufferIndex = 0;

    // For every boundary patch of the interface
    for (uint j = 0; j < patchIDs_.size(); j++)
    {
        int patchID = patchIDs_.at(j);

        // Get the pressureRgh gradient boundary patch
        const scalarField gradientPatch((pRgh_->boundaryFieldRef()[patchID]).snGrad());

        // For every cell of the patch
        forAll(gradientPatch, i)
        {
            // Copy the pressureRgh gradient into the buffer
            buffer[bufferIndex++] = -gradientPatch[i];
        }
    }
}

void preciceAdapter::FF::PressureRghGradient::read(double* buffer, const unsigned int dim)
{
    int bufferIndex = 0;

    // For every boundary patch of the interface
    for (uint j = 0; j < patchIDs_.size(); j++)
    {
        int patchID = patchIDs_.at(j);

        // Get the pressureRgh gradient boundary patch
        scalarField& gradientPatch = refCast<fixedGradientFvPatchScalarField>(
                pRgh_->boundaryFieldRef()[patchID]).gradient();

        // For every cell of the patch
        forAll(gradientPatch, i)
        {
            // Set the pressureRgh gradient as the buffer value
            gradientPatch[i] = buffer[bufferIndex++];
        }
    }
}

bool preciceAdapter::FF::PressureRghGradient::isLocationTypeSupported(const bool meshConnectivity) const
{
    return (this->locationType_ == LocationType::faceCenters);
}

std::string preciceAdapter::FF::PressureRghGradient::getDataName() const
{
    return "PressureGradient";
}
