#ifndef robinBC2FvPatchScalarField_H
#define robinBC2FvPatchScalarField_H

#include "fvPatchFields.H"
#include "fvPatchField.H"
#include "fieldTypes.H"

namespace Foam
{
    class robinBC2FvPatchScalarField : public fvPatchScalarField
    {
        // Private data
        scalar lowerAngle_;
        scalar upperAngle_;
        scalar fixedValue_;

    public:
        //- Runtime type information
        TypeName("robinBC2");

        // Constructors

        //- Construct from patch and internal field
        robinBC2FvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        robinBC2FvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given robinBC2FvPatchScalarField
        //  onto a new patch
        robinBC2FvPatchScalarField
        (
            const robinBC2FvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        robinBC2FvPatchScalarField
        (
            const robinBC2FvPatchScalarField&
        );

        // Member functions

        //- Evaluate the boundary condition
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
    };
}

#endif

