/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2013 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
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

Description:
    This boundary condition employs the Kistler dynamic contact angle model (taken from Dirk Drunding and 
    Anja Lippert dissertation)

    // theta_d = f_Hoff (Ca + f^(-1)_Hoff (\theta0_))
    // f_Hoff (x) := arccos (1- 2*tanh(5.16(x / (1 + 1.31*x^(0.99)))^0.706))

Developed by:
    Muhammad Hassan Asghar
    Mathematical Modeling and Analysis
    TU  Darmstadt
\*---------------------------------------------------------------------------*/

#include "dynamicAlphaContactAngleKistlerFvPatchScalarField.H"

#include "addToRunTimeSelectionTable.H"

#include "fvPatchFieldMapper.H"

#include "volMesh.H"

#include "volFields.H"

#include "mathematicalConstants.H"

#include <cmath>

#include <functional>

#include "interpolationCellPoint.H"

// * * * * * * * * * * * Private Member Functions * * * * * * * * * * * * * * //

bool Foam::dynamicAlphaContactAngleKistlerFvPatchScalarField::
hasContactLine
(
    label faceI
) const
{
        // Look up PLIC normals and positions. 
    const auto& db = this->db(); 

    const auto normalsName = IOobject::groupName
    (
        "interfaceNormal", 
        this->internalField().group()
    );
    const auto centresName = IOobject::groupName
    (
        "interfaceCentre", 
        this->internalField().group()
    );

    bool hasNormals = db.found(normalsName);
    if (!hasNormals)
    {
        // This BC is updated before interface reconstruction.
        // Do nothing if PLIC fields are not available in the registry. 
        return false;
    }

    bool hasCentres = db.found(centresName);
    if (!hasCentres)
    {
        // This BC is updated before interface reconstruction.
        // Do nothing if PLIC fields are not available in the registry. 
        return false;
    }

    const volVectorField& interfaceNormal = 
        db.lookupObject<volVectorField>(normalsName);

    const volVectorField& interfaceCentre = 
        db.lookupObject<volVectorField>(centresName);

    // Get patch fields for interface normals and centers
    const fvPatch& patch = this->patch();
    const label patchIndex = patch.index();
    const auto& pInterfaceNormals = interfaceNormal.boundaryField()[patchIndex];
    const auto& pInterfaceCentres = interfaceCentre.boundaryField()[patchIndex];

    // Get patch internal fields of normals and centers
    const auto pInternalNormalsTmp = pInterfaceNormals.patchInternalField();
    const auto& pInternalNormals = pInternalNormalsTmp.cref(); 
    const auto pInternalCentresTmp = pInterfaceCentres.patchInternalField();
    const auto& pInternalCentres = pInternalCentresTmp.cref(); 

    const vector& cellInterfaceNormal = pInternalNormals[faceI];
    const vector& cellInterfaceCentre = pInternalCentres[faceI];

    const auto& mesh = interfaceNormal.mesh();
    const auto& meshPoints = mesh.points();
    const auto& meshFaces = mesh.faces();
    const auto& thisFace = meshFaces[patch.start() + faceI];

    // Get face points. 
    for(auto pointI = 0; pointI < (thisFace.size() - 1); ++pointI)
    {
        // Compute the signed distance of the first point.
        const point& firstFacePoint = meshPoints[thisFace[pointI]];
        const scalar firstDist = (firstFacePoint - cellInterfaceCentre) & 
            cellInterfaceNormal;

        // Compute the signed distance of the second point.
        const point& secondFacePoint = meshPoints[thisFace[pointI + 1]];
        const scalar secondDist = (secondFacePoint - cellInterfaceCentre) & 
            cellInterfaceNormal;

        if (firstDist * secondDist < 0)
        {
            return true;
        }
    }

    return false;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dynamicAlphaContactAngleKistlerFvPatchScalarField::
dynamicAlphaContactAngleKistlerFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    alphaContactAngleTwoPhaseFvPatchScalarField(p, iF),
    theta0_(0.0),
    contactLineAngle_
    (
        IOobject
        (
            "clangle", 
            this->patch().boundaryMesh().mesh().time().timeName(),  
            this->patch().boundaryMesh().mesh(), 
            IOobject::NO_READ, 
            IOobject::AUTO_WRITE
        ),
        this->patch().boundaryMesh().mesh(), 
        dimensionedScalar("clangle", dimless, 0)
    )
{}


Foam::dynamicAlphaContactAngleKistlerFvPatchScalarField::
dynamicAlphaContactAngleKistlerFvPatchScalarField
(
    const dynamicAlphaContactAngleKistlerFvPatchScalarField& gcpsf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    alphaContactAngleTwoPhaseFvPatchScalarField(gcpsf, p, iF, mapper),
    theta0_(gcpsf.theta0_),
    contactLineAngle_
    (
        IOobject
        (
            "clangle", 
            this->patch().boundaryMesh().mesh().time().timeName(),  
            this->patch().boundaryMesh().mesh(), 
            IOobject::NO_READ, 
            IOobject::AUTO_WRITE
        ),
        this->patch().boundaryMesh().mesh(), 
        dimensionedScalar("clangle", dimless, 0)
    )
{}


Foam::dynamicAlphaContactAngleKistlerFvPatchScalarField::
dynamicAlphaContactAngleKistlerFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    alphaContactAngleTwoPhaseFvPatchScalarField(p, iF, dict),
    theta0_(dict.get<scalar>("theta0")),
    contactLineAngle_
    (
        IOobject
        (
            "clangle", 
            this->patch().boundaryMesh().mesh().time().timeName(),  
            this->patch().boundaryMesh().mesh(), 
            IOobject::NO_READ, 
            IOobject::AUTO_WRITE
        ),
        this->patch().boundaryMesh().mesh(), 
        dimensionedScalar("clangle", dimless, 0)
    ) 
{
    evaluate();
}


Foam::dynamicAlphaContactAngleKistlerFvPatchScalarField::
dynamicAlphaContactAngleKistlerFvPatchScalarField
(
    const dynamicAlphaContactAngleKistlerFvPatchScalarField& gcpsf
)
:
    alphaContactAngleTwoPhaseFvPatchScalarField(gcpsf),
    theta0_(gcpsf.theta0_),
    contactLineAngle_
    (
        IOobject
        (
            "clangle", 
            this->patch().boundaryMesh().mesh().time().timeName(),  
            this->patch().boundaryMesh().mesh(), 
            IOobject::NO_READ, 
            IOobject::AUTO_WRITE
        ),
        this->patch().boundaryMesh().mesh(), 
        dimensionedScalar("clangle", dimless, 0)
    )
{}


Foam::dynamicAlphaContactAngleKistlerFvPatchScalarField::
dynamicAlphaContactAngleKistlerFvPatchScalarField
(
    const dynamicAlphaContactAngleKistlerFvPatchScalarField& gcpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    alphaContactAngleTwoPhaseFvPatchScalarField(gcpsf, iF),
    theta0_(gcpsf.theta0_),
    contactLineAngle_
    (
        IOobject
        (
            "clangle", 
            this->patch().boundaryMesh().mesh().time().timeName(),  
            this->patch().boundaryMesh().mesh(), 
            IOobject::NO_READ, 
            IOobject::AUTO_WRITE
        ),
        this->patch().boundaryMesh().mesh(), 
        dimensionedScalar("clangle", dimless, 0)
    ) 
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField>
Foam::dynamicAlphaContactAngleKistlerFvPatchScalarField::theta
(
    const fvPatchVectorField& Up,
    const fvsPatchVectorField& nHat
) const
{
    //Patch Index
    const auto& patch = this->patch();
    const label patchIndex =patch.index();

    //patch face normal vectors
    const vectorField nf(patch.nf());

    const dictionary& transportProperties =
    this->db().objectRegistry::lookupObject<IOdictionary>
    (
        "transportProperties"
    );

    // Choose the larger dynamic viscosity from available two phases.
    word phase1Name (wordList(transportProperties.lookup("phases"))[0]);
    word phase2Name (wordList(transportProperties.lookup("phases"))[1]);

    // Get constant phase-specific densities and kinematic viscosities.
    dimensionedScalar rho1(transportProperties.subDict(phase1Name).get<dimensionedScalar>("rho"));
    dimensionedScalar nu1c(transportProperties.subDict(phase1Name).get<dimensionedScalar>("nu"));

    dimensionedScalar rho2(transportProperties.subDict(phase2Name).get<dimensionedScalar>("rho"));
    dimensionedScalar nu2c(transportProperties.subDict(phase2Name).get<dimensionedScalar>("nu"));

    word nuName; 
    dimensionedScalar rho(rho1); 
    // If the dynamic viscosity of phase1 is larger
    if (rho1*nu1c > rho2*nu2c)
    {
        nuName = "nu1";
        rho = rho1;
    }
    else
    {
        nuName = "nu2";
        rho = rho2;
    }
    const volScalarField& nu =
        this->db().objectRegistry::lookupObject<volScalarField>(nuName);

    // Fetch the wall kinematic viscosity of phase1 
    const label patchi = this->patch().index();
    const fvPatchScalarField& nup = nu.boundaryField()[patchi];
    scalarField muwall (nup*rho.value());

    // Wall Capillary number
    dimensionedScalar sigmap(transportProperties.get<dimensionedScalar>("sigma"));

    word fieldName = "alpha." + phase1Name;
    const volScalarField& alpha1_ =
        nu.mesh().lookupObject<volScalarField>(fieldName);

    const volScalarField::Boundary& abf = alpha1_.boundaryField();   
    // Lookup the desired alpha values on the patch
    const fvPatchField<scalar>& alphaPatchField = abf[patchIndex];
    
    //Initialization without velocity which is updated later in the code
    //Capillary number fields
    scalarField Ca(muwall / sigmap.value());

    // 1. Get the PLIC centre and normals
    // Look up PLIC normals and positions. 
    const auto& db = this->db(); 

    const auto normalsName = IOobject::groupName
    (
        "interfaceNormal", 
        this->internalField().group()
    );
    const auto centresName = IOobject::groupName
    (
        "interfaceCentre", 
        this->internalField().group()
    );

    bool hasNormals = db.found(normalsName);
    bool hasCentres = db.found(centresName);
            
    // Cell-centred velocity field
    Foam::volVectorField& U  = const_cast<Foam::volVectorField&> (this -> db().objectRegistry::lookupObject < volVectorField > ("U"));

    // velocity field at the plic centre
    // is calculated using inverse distance interpolation scheme using U and x_c
    //This field is then used to calculate the Capillary number as in isoAdvection.C
    
    // Uplic is used for advection: interpolated velocity at the PLIC centre
    Foam::volVectorField Uplic (U*scalar(0));

    //If there are no PLIC elements i.e, no wall velocity making Capillary number equals zero
    // Else calculate Ca using interpolated cell centred velocity to PLIC centre (using inverse distance weightings)
    if (!hasNormals || !hasCentres )
    {
        Info << "not Interpolating the PLIC" << nl;
        // This BC is updated before interface reconstruction.
        // Do nothing if PLIC fields are not available in the registry. 

        // Calculated the component of the velocity parallel to the wall
        vectorField Uwall(Up.patchInternalField() - Up);
        Uwall -= (nf & Uwall) * nf;

        // Find the direction of the interface parallel to the wall
        vectorField nWall(nHat - (nf & nHat) * nf);

        // Normalise nWall
        nWall /= (mag(nWall) + SMALL);

        // Calculate Uwall resolved normal to the interface parallel to
        // the interface
        scalarField uwall(-nWall & Uwall);

        //Capillary number fields
        Ca = muwall * uwall / sigmap.value();
    }
    else
    {
        Info << "Interpolating the PLIC" << nl;
        const volVectorField& interfaceCentre = 
            db.lookupObject<volVectorField>(centresName);

        // Create object for interpolating velocity to isoface centres
        interpolationCellPoint<vector> UInterp(U);

        //Interpolate the cell-centred velocity to the PLIC-centre
        forAll(interfaceCentre, cellI)
        {
            Uplic[cellI] = UInterp.interpolate(interfaceCentre[cellI], cellI);
        }
        //Get patch fields for PLIC velocity
        const auto& pUplic = Uplic.boundaryField()[patchIndex];

        // Calculated the component of the velocity parallel to the wall
        vectorField Uwall(pUplic.patchInternalField() - pUplic);

        Uwall -= (nf & Uwall) * nf;

        // Find the direction of the interface parallel to the wall
        vectorField nWall(nHat - (nf & nHat) * nf);

        // Normalise nWall
        nWall /= (mag(nWall) + SMALL);

        // Calculate Uwall resolved normal to the interface parallel to
        // the interface
        scalarField uwall(-nWall & Uwall);

        //Capillary number fields
        Ca = muwall * uwall / sigmap.value();
    }

    // Compute the contact angles at the wall.
    tmp<scalarField> thetafTmp = Foam::radToDeg(Foam::acos(nHat & nf));
    scalarField& thetaf = thetafTmp.ref();   

    // Visualization.
    auto& clangleBoundaryField = contactLineAngle_.boundaryFieldRef();
    auto& clanglePatchField = clangleBoundaryField[patchIndex];
    clanglePatchField *= 0.0;
    contactLineAngle_*=0.0;
        
    // Visualization of the contact angle
    const fvMesh& mesh = nu.mesh(); 
    const auto& faceOwner = mesh.faceOwner();   

    
    forAll(thetaf, faceI)
    {   
        if (hasContactLine(faceI))
        {
            if (hasNormals)
            {
                const volVectorField& interfaceCentre = 
                    db.lookupObject<volVectorField>(centresName);

                const volVectorField& interfaceNormal = 
                    db.lookupObject<volVectorField>(normalsName);

                // Visualize current PLIC contact angle as a cell-centered value.
                scalar cellI = faceOwner[patch.start() + faceI];
                scalar contactAngle = acos((interfaceNormal[cellI]) & nf[faceI] / (mag (interfaceNormal[cellI]) * mag(nf[faceI])) ) * 180.0 / M_PI;
                contactLineAngle_[faceOwner[patch.start() + faceI]] = contactAngle;
                
            }
            else
            {
                    // Visualize current PLIC contact angle as a cell-centered value.
                    contactLineAngle_[faceOwner[patch.start() + faceI]] = thetaf[faceI];
            }

            thetaf[faceI] = min(
                                180 / constant::mathematical::pi * kistler_model(theta0_, mag(Ca[faceI])),
                                scalar(180)
                                ); 
            clanglePatchField[faceI] = thetaf[faceI];

        }
    }
    return thetafTmp;
}
// Newton's method for root finding
double Foam::dynamicAlphaContactAngleKistlerFvPatchScalarField::newton_method(const std::function<double(double)>& func, const std::function<double(double)>& deriv, scalar initial_guess) const
{
    scalar x = initial_guess;
    scalar h = func(x) / deriv(x);
    while (std::abs(h) >= 1e-6) {
        h = func(x) / deriv(x);
        x = x - h;
    }
    return x;
}

//Inverse Hoffmann function
double Foam::dynamicAlphaContactAngleKistlerFvPatchScalarField::hoffman_function_inv(scalar theta0_) const
{
    scalar theta0_rad = theta0_* M_PI / 180.0;
    scalar right_side = std::pow((std::atanh((1 - std::cos(theta0_rad)) / 2) / 5.16), 1 / 0.706);

    auto equation = [right_side](double y) -> double {
        return y - (1.31 * right_side * std::pow(y, 0.99)) - right_side;
    };

    auto derivative = [right_side](double y) -> double {
        return 1 - (1.31 * 0.99 * right_side * std::pow(y, -0.01));
    };

    scalar initial_guess = 0.1;
    double reValue =newton_method(equation, derivative, initial_guess);
    return reValue;
}


// Kistler model
double Foam::dynamicAlphaContactAngleKistlerFvPatchScalarField::kistler_model(scalar theta0_, scalar Ca)  const 
{
    double invHoff = hoffman_function_inv(theta0_);
    double x = invHoff + mag(Ca);
    double theta_d = std::acos(1 - 2 * std::tanh(5.16 * std::pow(x / (1 + 1.31*std::pow(x, 0.99)), 0.706)));
    return theta_d;
}

void Foam::dynamicAlphaContactAngleKistlerFvPatchScalarField::write(Ostream& os) const
{
    alphaContactAngleTwoPhaseFvPatchScalarField::write(os);
    os.writeKeyword("theta0") << theta0_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        dynamicAlphaContactAngleKistlerFvPatchScalarField
    );
}


// ************************************************************************* //
