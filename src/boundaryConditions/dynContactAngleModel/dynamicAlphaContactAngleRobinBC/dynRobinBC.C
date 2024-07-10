#include "dynRobinBC.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "dictionary.H"
#include "alphaContactAngleTwoPhaseFvPatchScalarField.H"
#include "fvMesh.H"
#include "mathematicalConstants.H"
#include <cmath>
#include <functional>

namespace Foam
{

dynRobinBC::dynRobinBC
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    thetaA_(0.0),
    thetaR_(0.0),
    theta0_(0.0),
    initialPatchField_(p.size())
{
    // Store the initial patch volume fraction values
    forAll(*this, i)
    {
        initialPatchField_[i] = (*this)[i];
    }
}

dynRobinBC::dynRobinBC
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict),
    thetaA_(readScalar(dict.lookup("thetaA"))),
    thetaR_(readScalar(dict.lookup("thetaR"))),
    theta0_(readScalar(dict.lookup("theta0"))),
    initialPatchField_(p.size())
{
    // Store the initial patch volume fraction values
    forAll(*this, i)
    {
        initialPatchField_[i] = (*this)[i];
    }
}

dynRobinBC::dynRobinBC
(
    const dynRobinBC& bc,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(bc, p, iF, mapper),
    thetaA_(bc.thetaA_),
    thetaR_(bc.thetaR_),
    theta0_(bc.theta0_),
    initialPatchField_(bc.initialPatchField_)
{}

dynRobinBC::dynRobinBC(const dynRobinBC& bc)
:
    fixedValueFvPatchScalarField(bc),
    thetaA_(bc.thetaA_),
    thetaR_(bc.thetaR_),
    theta0_(bc.theta0_),
    initialPatchField_(bc.initialPatchField_)
{}

void dynRobinBC::autoMap(const fvPatchFieldMapper& m)
{
    fixedValueFvPatchScalarField::autoMap(m);
}

void dynRobinBC::rmap(const fvPatchScalarField& ptf, const labelList& addr)
{
    fixedValueFvPatchScalarField::rmap(ptf, addr);
}

void dynRobinBC::updateCoeffs()
{
    if (updated())
    {
        return;
    }
    

    const scalarField& patchField = *this;
    scalarField theta = patchField;

    const auto& mesh = this->patch().boundaryMesh().mesh();

    //patch face normal vectors
    const vectorField nf(patch().nf());

    if (mesh.time().timeIndex() <=1 )
    {
        // Inside hysteresis region, apply fixedValue BC
        forAll(patchField, i)
        {
            operator==(initialPatchField_[i]);
        }
    }
    else
    {

        const volVectorField& interfaceNormals_ = 
            this->db().lookupObject<volVectorField>("interfaceNormal");

        forAll(patchField, i)
        {
            scalar localTheta = acos(interfaceNormals_[patch().faceCells()[i]] & nf[i]) * 180 / M_PI;
            if (localTheta > thetaR_ && localTheta < thetaA_)
            {
                // Inside hysteresis region, apply fixedValue BC
                operator==(initialPatchField_[i]);
            }
            else
            {
                // Outside hysteresis region, apply Neumann BC
                operator==(0.0);
            }
        }
    }

    fixedValueFvPatchScalarField::updateCoeffs();
}

tmp<Foam::scalarField> dynRobinBC::theta
(
    const fvPatchVectorField& Up,
    const fvsPatchVectorField& nHat
) const
{
    scalarField localTheta(patch().size());

    //patch face normal vectors
    const vectorField nf(patch().nf());

    const auto& mesh = this->patch().boundaryMesh().mesh();

    if (mesh.time().timeIndex() <=1 )
    {
        return Foam::radToDeg(Foam::acos(nHat & nf));
    }
    else
    {
        const volVectorField& interfaceNormals_ = 
            this->db().lookupObject<volVectorField>("interfaceNormal");

        forAll(localTheta, i)
        {
            localTheta[i] = acos(interfaceNormals_[patch().faceCells()[i]] & nf[i]) * 180 /M_PI;
        }

        return tmp<scalarField>(new scalarField(localTheta));
    }

   
}

void dynRobinBC::write(Ostream& os) const
{
    fixedValueFvPatchScalarField::write(os);
    os.writeEntry("thetaA", thetaA_);
    os.writeEntry("thetaR", thetaR_);
    os.writeEntry("theta0", theta0_);
}

}

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        dynRobinBC
    );
}
