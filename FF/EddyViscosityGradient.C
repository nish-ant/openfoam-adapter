#include "EddyViscosityGradient.H"

using namespace Foam;

preciceAdapter::FF::EddyViscosityGradient::EddyViscosityGradient(
    const Foam::fvMesh& mesh,
    const std::string nameNut)
: nut_(
    const_cast<volScalarField*>(
        &mesh.lookupObject<volScalarField>(nameNut)))
{
    dataType_ = scalar;
}

void preciceAdapter::FF::EddyViscosityGradient::write(double* buffer, bool meshConnectivity, const unsigned int dim)
{
    int bufferIndex = 0;

    // For every boundary patch of the interface
    for (uint j = 0; j < patchIDs_.size(); j++)
    {
        int patchID = patchIDs_.at(j);

        // Get the eddy viscosity gradient boundary patch
        const scalarField gradientPatch((nut_->boundaryFieldRef()[patchID])
                                            .snGrad());

        // For every cell of the patch
        forAll(gradientPatch, i)
        {
            // Copy the eddy viscosity gradient into the buffer
            buffer[bufferIndex++] = -gradientPatch[i];
        }
    }
}

void preciceAdapter::FF::EddyViscosityGradient::read(double* buffer, const unsigned int dim)
{
    int bufferIndex = 0;

    // For every boundary patch of the interface
    for (uint j = 0; j < patchIDs_.size(); j++)
    {
        int patchID = patchIDs_.at(j);

        // Get the eddy viscosity gradient boundary patch
        scalarField& gradientPatch = refCast<fixedGradientFvPatchScalarField>(
                nut_->boundaryFieldRef()[patchID]).gradient();

        // For every cell of the patch
        forAll(gradientPatch, i)
        {
            // Set the eddy viscosity gradient as the buffer value
            gradientPatch[i] = buffer[bufferIndex++];
        }
    }
}

bool preciceAdapter::FF::EddyViscosityGradient::isLocationTypeSupported(const bool meshConnectivity) const
{
    return (this->locationType_ == LocationType::faceCenters);
}

std::string preciceAdapter::FF::EddyViscosityGradient::getDataName() const
{
    return "EddyViscosityGradient";
}
