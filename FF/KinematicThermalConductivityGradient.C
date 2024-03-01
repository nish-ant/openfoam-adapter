#include "KinematicThermalConductivityGradient.H"

using namespace Foam;

preciceAdapter::FF::KinematicThermalConductivityGradient::KinematicThermalConductivityGradient(
    const Foam::fvMesh& mesh,
    const std::string nameKappat)
: kappat_(
    const_cast<volScalarField*>(
        &mesh.lookupObject<volScalarField>(nameKappat)))
{
    dataType_ = scalar;
}

void preciceAdapter::FF::KinematicThermalConductivityGradient::write(double* buffer, bool meshConnectivity, const unsigned int dim)
{
    int bufferIndex = 0;

    // For every boundary patch of the interface
    for (uint j = 0; j < patchIDs_.size(); j++)
    {
        int patchID = patchIDs_.at(j);

        // Get the kinematic thermal conductivity gradient boundary patch
        const scalarField gradientPatch((kappat_->boundaryFieldRef()[patchID]).snGrad());

        // For every cell of the patch
        forAll(gradientPatch, i)
        {
            // Copy the kinematic thermal conductivity gradient into the buffer
            buffer[bufferIndex++] = -gradientPatch[i];
        }
    }
}

void preciceAdapter::FF::KinematicThermalConductivityGradient::read(double* buffer, const unsigned int dim)
{
    int bufferIndex = 0;

    // For every boundary patch of the interface
    for (uint j = 0; j < patchIDs_.size(); j++)
    {
        int patchID = patchIDs_.at(j);

        // Get the kinematic thermal conductivity gradient boundary patch
        scalarField& gradientPatch = refCast<fixedGradientFvPatchScalarField>(
                kappat_->boundaryFieldRef()[patchID]).gradient();

        // For every cell of the patch
        forAll(gradientPatch, i)
        {
            // Set the kinematic thermal conductivity gradient as the buffer value
            gradientPatch[i] = buffer[bufferIndex++];
        }
    }
}

bool preciceAdapter::FF::KinematicThermalConductivityGradient::isLocationTypeSupported(const bool meshConnectivity) const
{
    return (this->locationType_ == LocationType::faceCenters);
}

std::string preciceAdapter::FF::KinematicThermalConductivityGradient::getDataName() const
{
    return "KinematicThermalConductivityGradient";
}
