/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017-2021 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "navierSlipHysteresis.H"
#include "symmTransformField.H"
#include "fvPatchFieldMapper.H"

#include "volMesh.H"

#include "volFields.H"

#include "mathematicalConstants.H"

#include <cmath>

#include <functional>

#include "interpolationCellPoint.H"

#include "fvcDdt.H"
#include "fvcDiv.H"
#include "fvcLaplacian.H"
#include "fvcReconstruct.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::navierSlipHysteresis<Type>::navierSlipHysteresis
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    parent_bctype(p, iF),
    refValue_(p.size(), Zero),
    slipLength_(p.size(), 0.0001),
    gamma_(0.0),
    writeValue_(false)
{}


template<class Type>
Foam::navierSlipHysteresis<Type>::navierSlipHysteresis
(
    const navierSlipHysteresis<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    parent_bctype(ptf, p, iF, mapper),
    refValue_(ptf.refValue_, mapper),
    slipLength_(ptf.slipLength_, mapper),
    gamma_(ptf.gamma_),
    writeValue_(ptf.writeValue_)
{}


template<class Type>
Foam::navierSlipHysteresis<Type>::navierSlipHysteresis
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    parent_bctype(p, iF),
    refValue_(p.size(), Zero),
    slipLength_("slipLength", dict, p.size()),
    gamma_(dict.get<scalar>("gamma")),
    writeValue_(dict.getOrDefault("writeValue", false))
{
    fvPatchFieldBase::readDict(dict);

    // Backwards compatibility - leave refValue as zero unless specified
    refValue_.assign("refValue", dict, p.size(), IOobjectOption::LAZY_READ);

    evaluate();
}


template<class Type>
Foam::navierSlipHysteresis<Type>::navierSlipHysteresis
(
    const navierSlipHysteresis<Type>& ptf
)
:
    parent_bctype(ptf),
    refValue_(ptf.refValue_),
    slipLength_(ptf.slipLength_),
    gamma_(ptf.gamma_),
    writeValue_(ptf.writeValue_)
{}


template<class Type>
Foam::navierSlipHysteresis<Type>::navierSlipHysteresis
(
    const navierSlipHysteresis<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    parent_bctype(ptf, iF),
    refValue_(ptf.refValue_),
    slipLength_(ptf.slipLength_),
    gamma_(ptf.gamma_),
    writeValue_(ptf.writeValue_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

#include "fvPatchFields.H"
#include "volFields.H"
#include "dimensionedScalar.H"
#include <cmath>

// Function to compute f(x) for viscous stress
template<class Type>
Foam::scalarField Foam::navierSlipHysteresis<Type>::fStress(
    const vectorField& tangentialStress, // stress values for the patch
    const dimensionedScalar gamma, // limit for slip
    const scalar f_max, // maximum value of the function = 1.0
    const dimensionedScalar x_s, // the value of Stress where the function f=f_max; x_s = gamma +\delta
    const scalar n // exponent value to control the saturation position, must be b/w 1 and 2 to get monotonicity.
)
{
    // Initialize a scalarField for the function values
    scalarField fValues(tangentialStress.size(), 0.0);

    // Loop through the patchStressField to compute f(x) for each value
    forAll(tangentialStress, i)
    {
        scalar stress = mag(tangentialStress[i]); // Current stress value

        if (stress < gamma.value())
        {
            // If stress is below the threshold, f(x) = 0
            fValues[i] = 0.0;
        }
        else if (stress < x_s.value())
        {
            // If stress is between gamma and x_s, compute f(x) using the formula
            scalar num = pow(stress - gamma.value(), n);
            scalar denom = num + pow(x_s.value() - stress, n);
            fValues[i] = f_max * num / denom;
        }
        else
        {
            // If stress is above x_s, f(x) = f_max
            fValues[i] = f_max;
        }
    }

    return fValues; // Return the computed scalarField
}

template<class Type>
Foam::scalarField Foam::navierSlipHysteresis<Type>::computeValueFraction() const
{
    const fvPatch& patch = this->patch();
    const label patchIndex = patch.index();

    // Patch face normal vectors
    const vectorField nf(patch.nf());

    
     // // Look up PLIC normals and positions. 
    const objectRegistry& db = this->db(); // Ensure db is defined correctly

    //Get kinematic viscosity
    const volScalarField& nu =
        this->db().objectRegistry::lookupObject<volScalarField>("nu");
    const auto& mesh = nu.mesh();

    // Get the velocity field and compute the gradient
    const volVectorField& U_= mesh.lookupObject<volVectorField>("U");
    volTensorField gradU = fvc::grad(U_, "pointCellsLeastSquares");

    // Dynamic viscosity
    dimensionedScalar mu(
        "mu",
        dimensionSet(1, -1, -1, 0, 0, 0, 0), 
        1e-3 
    );

    // Compute viscous stress
    volSymmTensorField viscousStress = mu * (gradU + gradU.T());
    vectorField normalStress = viscousStress.boundaryField()[patchIndex] & nf;

    // Tangential projection tensor
    tensorField Ptan = I - sqr(nf);
    vectorField tangentialStress = transform(Ptan, normalStress);

    // Compute f(|S|)
    dimensionedScalar gamma(
        "gamma",
        dimensionSet(1, -1, -2, 0, 0, 0, 0),
        1e-3  // Slip threshold
    );
    dimensionedScalar x_s(
        "x_s",
        gamma.dimensions(),
        1.1 * gamma.value()  // Saturation value for f(S)
    );

    scalarField fValues = fStress(tangentialStress, gamma, 1.0, x_s, 2);

    // Compute valueFraction
    scalar d = 0.0001;  // Cell size
    scalarField valueFraction(patch.size());

    forAll(valueFraction, i)
    {
        valueFraction[i] = d / (d + (2 * slipLength[i] * fValues[i]));
    }

    return valueFraction;
}


template<class Type>
void Foam::navierSlipHysteresis<Type>::autoMap
(
    const fvPatchFieldMapper& m
)
{
    parent_bctype::autoMap(m);
    refValue_.autoMap(m);
    slipLength_.autoMap(m);
}


template<class Type>
void Foam::navierSlipHysteresis<Type>::rmap
(
    const fvPatchField<Type>& ptf,
    const labelList& addr
)
{
    parent_bctype::rmap(ptf, addr);

    const auto& dmptf =
        refCast<const navierSlipHysteresis<Type>>(ptf);

    refValue_.rmap(dmptf.refValue_, addr);
    slipLength_.rmap(dmptf.slipLength_, addr);
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::navierSlipHysteresis<Type>::snGrad() const
{
    tmp<vectorField> nHat = this->patch().nf();
    const Field<Type> pif(this->patchInternalField());

    // Compute the modified valueFraction
    Foam::scalarField valueFraction = computeValueFraction();

    return
    (
        lerp
        (
            transform(I - sqr(nHat), pif),
            refValue_,
            valueFraction
        ) - pif
    )*this->patch().deltaCoeffs();
}

// Boundary condition code
template<class Type>
void Foam::navierSlipHysteresis<Type>::evaluate
(
    const Pstream::commsTypes
)
{
    if (!this->updated())
    {
        // Updates the boundary condition coefficient
        this->updateCoeffs();
    }


    // const auto normalsName = 
    //     IOobject::groupName("interfaceNormal", this->internalField().group());
    // const auto centresName = 
    //     IOobject::groupName("interfaceCentre", this->internalField().group());

    // bool hasNormals = db.found(normalsName);
    // bool hasCentres = db.found(centresName);

    // const volVectorField& interfaceNormal = 
    //     db.lookupObject<volVectorField>(normalsName);

    // const volVectorField& interfaceCentre = 
    //     db.lookupObject<volVectorField>(centresName);

    // // Get patch fields for interface normals and centers
    const fvPatch& patch = this->patch();
    const label patchIndex = patch.index();
    // const auto& pInterfaceNormals = interfaceNormal.boundaryField()[patchIndex];
    // const auto& pInterfaceCentres = interfaceCentre.boundaryField()[patchIndex];

    //patch face normal vectors
    const vectorField nf(patch.nf());
    
    // // Compute the contact angles at the wall.
    // tmp<scalarField> thetafTmp = Foam::radToDeg(Foam::acos(pInterfaceNormals & nf));
    // scalarField& thetaf = thetafTmp.ref(); 

    // scalar theta0 = 45;

    // // Fetch surface tension coefficient
    // const dictionary& transportProperties =
    // this->db().objectRegistry::lookupObject<IOdictionary>
    // (
    //     "transportProperties"
    // );
    // dimensionedScalar sigmap(transportProperties.get<dimensionedScalar>("sigma"));
    
    // // Compute Young's stress term
    // scalarField contactAngleDiff = cos(thetaf) - cos(theta0);  
    // dimensionedScalar epsilon
    // (
    //     "epsilon",  
    //     dimensionSet(1, 0, 0, 0, 0, 0, 0),  
    //     1e-9  // Value of the scalar 
    // );
    // // forAll(thetaf, faceI)
    // // {
    // //     youngsStress[faceI] = sigmap * contactAngleDiff *pInterfaceNormals[faceI] /epsilon;
    // // }
    // vectorField youngsStress = sigmap.value() * contactAngleDiff * pInterfaceNormals / epsilon.value();  


    //Calculate totalStress 
   // calculate the limit gamma of f(S)
    //calculate the funcation of stress
    // change he value fraction

    //totalStress = viscousStress + wettingStress
    //For single phase flow, totalStress = viscousStress

   // const auto& mesh = nu.mesh();
   // const volVectorField& U_= mesh.lookupObject<volVectorField>("U");

    //volTensorField gradU = fvc::grad(U_, "pointCellsLeastSquares"); 

    // Properties for water
    // dimensionedScalar mu
    // (
    //     "mu",  
    //     dimensionSet(1, -1, -1, 0, 0, 0, 0),  
    //     1e-3  // Value of the scalar 
    // );

    // volSymmTensorField viscousStress = mu * (gradU + gradU.T());
    // // Get viscous stress at the patch S.nf
    // vectorField normalStress = viscousStress.boundaryField()[patchIndex] & nf;

    // tensorField Ptan = I - sqr(nf);    // Tangential projection tensor

    // // Project the normal stress onto the tangential plane
    // vectorField tangentialStress = transform(Ptan, normalStress);

    // // Calculation of gamma; limit for slip
    // dimensionedScalar gamma
    // (
    //     "gamma",  
    //     dimensionSet(1, -1, -2, 0, 0, 0, 0),  
    //     1e-3  // Value of the scalar 
    // );

    // // stress values when the function saturates to max
    // Foam::dimensionedScalar x_s(
    //     "x_s",
    //     gamma.dimensions(),  // stress
    //     1.1 * gamma.value()    // Compute x_s = k * gamma
    // );

    // Foam::scalarField fValues = fStress(tangentialStress, gamma, 1.0, x_s, 2);

    // //Value fraction for each cell as it depends on the f(S)
    // Foam::scalarField valueFraction(this->patch().size());

    // Foam::scalar slipLength = 0.0001;
    // Foam::scalar d = 0.0001; //cell size

    // forAll(valueFraction, i)
    // {
    //     // Update the scaling factor 
    //     valueFraction[i] = d / (d + (2 * slipLength* fValues[i]));
    // }

    // Compute the scaling factor (valueFraction)
    scalarField valueFraction = computeValueFraction();

    Field<Type>::operator=
    (
        lerp
        (
            transform(I - sqr(this->patch().nf()), this->patchInternalField()),
            refValue_,
            valueFraction
        )
        //+ youngsStress
    );

    parent_bctype::evaluate();
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::navierSlipHysteresis<Type>::snGradTransformDiag() const
{
    tmp<vectorField> diag(cmptMag(this->patch().nf()));

    // Compute the modified valueFraction
    scalarField valueFraction = computeValueFraction();

    return
        valueFraction*pTraits<Type>::one
      + (1.0 - valueFraction)
       *transformFieldMask<Type>(pow<vector, pTraits<Type>::rank>(diag));
}


template<class Type>
void Foam::navierSlipHysteresis<Type>::write(Ostream& os) const
{
    this->parent_bctype::write(os);
    refValue_.writeEntry("refValue", os);
    slipLength_.writeEntry("slipLength", os);

    if (writeValue_)
    {
        os.writeEntry("writeValue", "true");
        fvPatchField<Type>::writeValueEntry(os);
    }
}


// ************************************************************************* //
