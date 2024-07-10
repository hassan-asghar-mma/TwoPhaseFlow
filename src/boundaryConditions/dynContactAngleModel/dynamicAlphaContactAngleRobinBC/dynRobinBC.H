#ifndef dynRobinBC_H
#define dynRobinBC_H

#include "fixedValueFvPatchFields.H"
#include "volFields.H"
#include "vectorField.H"
#include "alphaContactAngleTwoPhaseFvPatchScalarField.H"
#include "fvMesh.H"

namespace Foam
{

class dynRobinBC : public fixedValueFvPatchScalarField
{
private:
    scalar thetaA_; // Advancing contact angle
    scalar thetaR_; // Receding contact angle
    scalar theta0_; // Equilibrium contact angle
    scalarField initialPatchField_; // Initial patch volume fraction values


public:
    TypeName("dynRobinBC");

    dynRobinBC(
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&
    );

    dynRobinBC(
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&,
        const dictionary&
    );

    dynRobinBC(
        const dynRobinBC&,
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&,
        const fvPatchFieldMapper&
    );

    dynRobinBC(
        const dynRobinBC&
    );

    //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new dynRobinBC(*this)
            );
        }

        //- Construct as copy setting internal field reference
        dynRobinBC
        (
            const dynRobinBC&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new dynRobinBC(*this, iF)
            );
        }



    void autoMap(const fvPatchFieldMapper&) override;

    void rmap(const fvPatchScalarField&, const labelList&) override;

    void updateCoeffs() override;

    // Return the equilibrium contact-angle
    virtual tmp<scalarField> theta
    (
        const fvPatchVectorField& Up,
        const fvsPatchVectorField& nHat
    ) const;

    // Write
    virtual void write(Ostream&) const override;
};

}

#endif
