#include "TurbulentDissipationRateGradient.H"

using namespace Foam;

preciceAdapter::FF::TurbulentDissipationRateGradient::TurbulentDissipationRateGradient(
    const Foam::fvMesh& mesh,
    const std::string nameEpsilon)
: epsilon_(
    const_cast<volScalarField*>(
        &mesh.lookupObject<volScalarField>(nameEpsilon)))
{
    dataType_ = scalar;
}

void preciceAdapter::FF::TurbulentDissipationRateGradient::write(double* buffer, bool meshConnectivity, const unsigned int dim)
{
    int bufferIndex = 0;

    // For every boundary patch of the interface
    for (uint j = 0; j < patchIDs_.size(); j++)
    {
        int patchID = patchIDs_.at(j);

        // Get the Turbulent Dissipation Rate gradient boundary patch
        const scalarField gradientPatch((epsilon_->boundaryFieldRef()[patchID]).snGrad());

        // For every cell of the patch
        forAll(gradientPatch, i)
        {
            // Copy the Turbulent Dissipation Rate gradient into the buffer
            buffer[bufferIndex++] = -gradientPatch[i];
        }
    }
}

void preciceAdapter::FF::TurbulentDissipationRateGradient::read(double* buffer, const unsigned int dim)
{
    int bufferIndex = 0;

    // For every boundary patch of the interface
    for (uint j = 0; j < patchIDs_.size(); j++)
    {
        int patchID = patchIDs_.at(j);

        // Get the Turbulent Dissipation Rate gradient boundary patch
        scalarField& gradientPatch = refCast<fixedGradientFvPatchScalarField>(
                epsilon_->boundaryFieldRef()[patchID]).gradient();

        // For every cell of the patch
        forAll(gradientPatch, i)
        {
            // Set the Turbulent Dissipation Rate gradient as the buffer value
            gradientPatch[i] = buffer[bufferIndex++];
        }
    }
}

bool preciceAdapter::FF::TurbulentDissipationRateGradient::isLocationTypeSupported(const bool meshConnectivity) const
{
    return (this->locationType_ == LocationType::faceCenters);
}

std::string preciceAdapter::FF::TurbulentDissipationRateGradient::getDataName() const
{
    return "TurbulentDissipationRateGradient";
}
