#include "robinBC2FvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    // Runtime type information
    defineTypeNameAndDebug(robinBC2FvPatchScalarField, 0);

    addToRunTimeSelectionTable
    (
        fvPatchScalarField,
        robinBC2FvPatchScalarField,
        dictionary
    );

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    robinBC2FvPatchScalarField::robinBC2FvPatchScalarField
    (
        const fvPatch& p,
        const DimensionedField<scalar, volMesh>& iF
    )
    :
        fvPatchScalarField(p, iF),
        lowerAngle_(0),
        upperAngle_(0),
        fixedValue_(0)
    {}

    robinBC2FvPatchScalarField::robinBC2FvPatchScalarField
    (
        const fvPatch& p,
        const DimensionedField<scalar, volMesh>& iF,
        const dictionary& dict
    )
    :
        fvPatchScalarField(p, iF, dict)
    {
        lowerAngle_ = dict.lookupOrDefault<scalar>("lowerAngle", 0);
        upperAngle_ = dict.lookupOrDefault<scalar>("upperAngle", 90);
        fixedValue_ = dict.lookupOrDefault<scalar>("fixedValue", 0);
    }

    robinBC2FvPatchScalarField::robinBC2FvPatchScalarField
    (
        const robinBC2FvPatchScalarField& f,
        const fvPatch& p,
        const DimensionedField<scalar, volMesh>& iF,
        const fvPatchFieldMapper& m
    )
    :
        fvPatchScalarField(f, p, iF, m),
        lowerAngle_(f.lowerAngle_),
        upperAngle_(f.upperAngle_),
        fixedValue_(f.fixedValue_)
    {}

    robinBC2FvPatchScalarField::robinBC2FvPatchScalarField
    (
        const robinBC2FvPatchScalarField& f
    )
    :
        fvPatchScalarField(f),
        lowerAngle_(f.lowerAngle_),
        upperAngle_(f.upperAngle_),
        fixedValue_(f.fixedValue_)
    {}

    void robinBC2FvPatchScalarField::updateCoeffs()
    {
        if (updated())
        {
            return;
        }

        // Assume contactAngle is already available in the dictionary
        scalar contactAngle = 70;

        if (contactAngle >= lowerAngle_ && contactAngle <= upperAngle_)
        {
            // Fixed value boundary condition
            operator== (fixedValue_);
        }
        else
        {
            // Neumann boundary condition (zeroGradient)
            operator== (this->patchInternalField());
        }

        fvPatchScalarField::updateCoeffs();
    }

    void robinBC2FvPatchScalarField::write(Ostream& os) const
    {
        fvPatchScalarField::write(os);
        os.writeKeyword("lowerAngle") << lowerAngle_ << token::END_STATEMENT << nl;
        os.writeKeyword("upperAngle") << upperAngle_ << token::END_STATEMENT << nl;
        os.writeKeyword("fixedValue") << fixedValue_ << token::END_STATEMENT << nl;
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
}

