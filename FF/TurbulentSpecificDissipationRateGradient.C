#include "TurbulentSpecificDissipationRateGradient.H"

using namespace Foam;

preciceAdapter::FF::TurbulentSpecificDissipationRateGradient::TurbulentSpecificDissipationRateGradient(
    const Foam::fvMesh& mesh,
    const std::string nameOmega)
: omega_(
    const_cast<volScalarField*>(
        &mesh.lookupObject<volScalarField>(nameOmega)))
{
    dataType_ = scalar;
}

void preciceAdapter::FF::TurbulentSpecificDissipationRateGradient::write(double* buffer, bool meshConnectivity, const unsigned int dim)
{
    int bufferIndex = 0;

    // For every boundary patch of the interface
    for (uint j = 0; j < patchIDs_.size(); j++)
    {
        int patchID = patchIDs_.at(j);

        // Get the Turbulent Specific Dissipation Rate gradient boundary patch
        const scalarField gradientPatch((omega_->boundaryFieldRef()[patchID]).snGrad());

        // For every cell of the patch
        forAll(gradientPatch, i)
        {
            // Copy the Turbulent Specific Dissipation Rate gradient into the buffer
            buffer[bufferIndex++] = -gradientPatch[i];
        }
    }
}

void preciceAdapter::FF::TurbulentSpecificDissipationRateGradient::read(double* buffer, const unsigned int dim)
{
    int bufferIndex = 0;

    // For every boundary patch of the interface
    for (uint j = 0; j < patchIDs_.size(); j++)
    {
        int patchID = patchIDs_.at(j);

        // Get the Turbulent Specific Dissipation Rate gradient boundary patch
        scalarField& gradientPatch = refCast<fixedGradientFvPatchScalarField>(
                omega_->boundaryFieldRef()[patchID]).gradient();

        // For every cell of the patch
        forAll(gradientPatch, i)
        {
            // Set the Turbulent Specific Dissipation Rate gradient as the buffer value
            gradientPatch[i] = buffer[bufferIndex++];
        }
    }
}

bool preciceAdapter::FF::TurbulentSpecificDissipationRateGradient::isLocationTypeSupported(const bool meshConnectivity) const
{
    return (this->locationType_ == LocationType::faceCenters);
}

std::string preciceAdapter::FF::TurbulentSpecificDissipationRateGradient::getDataName() const
{
    return "TurbulentSpecificDissipationRateGradient";
}
