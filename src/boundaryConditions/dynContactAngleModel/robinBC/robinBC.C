/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
     \\/     M anipulation  |
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

#include "dimensionedScalarFwd.H"
#include "robinBC.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volMesh.H"
#include "fvCFD.H"
#include "volFields.H"
#include "mathematicalConstants.H"
#include "interpolationCellPoint.H"
#include <iostream>
#include <fstream>

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::robinBC::
robinBC
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    alphaContactAngleTwoPhaseFvPatchScalarField(p, iF),
    theta0_(0.0),
    ct_(0.0),
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


Foam::robinBC::
robinBC
(
    const robinBC& gcpsf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    alphaContactAngleTwoPhaseFvPatchScalarField(gcpsf, p, iF, mapper),
    theta0_(gcpsf.theta0_),
    ct_(gcpsf.ct_),
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


Foam::robinBC::
robinBC
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    alphaContactAngleTwoPhaseFvPatchScalarField(p, iF, dict),
    theta0_(readScalar(dict.lookup("theta0"))),
    ct_(readScalar(dict.lookup("ct"))),
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


Foam::robinBC::
robinBC
(
    const robinBC& gcpsf
)
:
    alphaContactAngleTwoPhaseFvPatchScalarField(gcpsf),
    theta0_(gcpsf.theta0_),
    ct_(gcpsf.ct_),
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


Foam::robinBC::
robinBC
(
    const robinBC& gcpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    alphaContactAngleTwoPhaseFvPatchScalarField(gcpsf, iF),
    theta0_(gcpsf.theta0_),
    ct_(gcpsf.ct_),
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

bool Foam::robinBC::
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

Foam::tmp<Foam::scalarField>
Foam::robinBC::theta
(
    const fvPatchVectorField& Up,
    const fvsPatchVectorField& nHat
) const
{

    const vectorField nf(patch().nf());

    // Calculated the component of the velocity parallel to the wall
    vectorField Uwall(Up.patchInternalField() - Up);
    Uwall -= (nf & Uwall)*nf;

    // Find the direction of the interface parallel to the wall
    vectorField nWall(nHat - (nf & nHat)*nf);

    // Normalise nWall
    nWall /= (mag(nWall) + SMALL);

    // Calculate Uwall resolved normal to the interface parallel to
    // the interface
    scalarField uwall(-nWall & Uwall);

    const auto& patch = this->patch();
    const label patchIndex = patch.index();

    const volScalarField& nu1 =
        this->db().objectRegistry::lookupObject<volScalarField>("nu1");

    const dictionary& transportProperties =
    this->db().objectRegistry::lookupObject<IOdictionary>
    (
        "transportProperties"
    );

    word phase1Name (wordList(transportProperties.lookup("phases"))[0]);
    word phase2Name (wordList(transportProperties.lookup("phases"))[1]);

    const auto& mesh_ = nu1.mesh();
    word fieldName = "alpha." + phase1Name;
    const volScalarField& alpha1_ = mesh_.lookupObject<volScalarField>(fieldName);
    Foam::volVectorField& U_  = const_cast<Foam::volVectorField&> (this -> db().objectRegistry::lookupObject < volVectorField > ("U"));

    dimensionedScalar rho1(transportProperties.subDict(phase1Name).get<dimensionedScalar>("rho"));
    dimensionedScalar sigmap(transportProperties.get<dimensionedScalar>("sigma"));
    const fvPatchScalarField&  nu1p = nu1.boundaryField()[patchIndex];

    scalarField Ca(nu1p*rho1.value()*uwall/sigmap.value());

    // thetaD = (ct*  Ca + theta0^3)^(1/3)

    //  Ca^1/3 in rad


//  steps:
// 1. Read the alpha_patch values and store them in the file -> done for the first time step only
// 2. Check ifthetaD is in hysteresis -> apply Dirichlet boundary condition -> reading the alpha values from the file
// 3. If outside proceed as normal    

    scalar thetaA_ = 110;
    scalar thetaR_ = 60;
    scalar theta0 = 90;
    // scalar to select between Dirichlet (1) or Neumann BC (0)
    scalar beta = 1;
    //Alpha field values at the patch i 
    auto& pAlphaField = const_cast<fvPatchScalarField&>(alpha1_.boundaryField()[patchIndex]);
    auto& pUField = const_cast<fvPatchVectorField&>(U_.boundaryField()[patchIndex]);

    // plicU is used for advection
    Foam::volVectorField Uplic (U_*scalar(0));
    
    if(mesh_.time().timeIndex()!=0)
    {
            // 1. Get the PLIC centre and normals
        // Look up PLIC normals and positions. 
        const auto& db = this->db(); 
        const auto centresName = IOobject::groupName
        (
            "interfaceCentre", 
            this->internalField().group()
        );

        // Create object for interpolating velocity to isoface centres
        interpolationCellPoint<vector> UInterp(U_);
        const volVectorField& interfaceCentre = db.lookupObject<volVectorField>(centresName);

        //Interpolate the cell-centred velocity to the PLIC-centre
        forAll(interfaceCentre, cellI)
        {
            Uplic[cellI] = UInterp.interpolate(interfaceCentre[cellI], cellI);
        }

    }
    //Get patch fields for PLIC velocity
    auto& pUplic = const_cast<fvPatchVectorField&>(Uplic.boundaryField()[patchIndex]);


    // Writing the patch alpha field into a file to read it again if inside hysteresis - Done only for the first time step
    // Open a file for writing
    if (mesh_.time().timeIndex() < 1)
    {
            std::ofstream outFile("patchAlphaField.csv");
           // std::ofstream UFile("patchUField.csv");
            // Check if the file is open
            if (outFile.is_open()) 
            {
                // Loop through the array and write each element to the file
                forAll(pAlphaField, i) {
                    outFile << pAlphaField[i] << std::endl; // Writing each element on a new line
                }
                // Close the file
                outFile.close();
            }
            // if(UFile.is_open())
            // {
            //     for (const auto& innerVec : pUField) {
            //         for (double num : innerVec) {
            //             UFile << num << " ";
            //         }
            //         UFile << "\n"; // New line for each inner vector
            //     }
            //     UFile.close();
            // }
    }

    std::vector <double> pAlphaFieldFile; //Values from t=0 read from the file 
    // std::vector<std::vector<double>> pUFieldFile; //Values from t=0 read from the file 

    std::ifstream inFile("patchAlphaField.csv");
    scalar value;
    if (inFile.is_open()) {
        while (inFile >> value) {
            pAlphaFieldFile.push_back(value);
        }
        inFile.close();
    }

    // std::ifstream inUFile("patchUField.csv");
    // if (inFile.is_open()) {
    //     std::string line;
    //     while (getline(inUFile, line)) {
    //         std::istringstream iss(line);
    //         std::vector<double> row;
    //         double value;

    //         while (iss >> value) {
    //             row.push_back(value);
    //         }
    //         pUFieldFile.push_back(row);
    //     }
    // }


    // Calculate cell-centered gradients at the patch i
    volVectorField gradAlpha = fvc::grad(alpha1_);
    const auto& pGradAlpha = gradAlpha.boundaryField()[patchIndex];
    
    //For calculation of distance between cell centre and face centre
    const volVectorField & C = mesh_.C(); // Cell center coordinates
    const surfaceVectorField & Cf = mesh_.Cf(); // Face center coordinates
    const auto & faceOwner = mesh_.faceOwner();

    // Compute the contact angles at the wall.
    tmp<scalarField> thetafTmp = Foam::radToDeg(Foam::acos(nHat & nf));
    scalarField& thetaf = thetafTmp.ref(); 
    //Info << "Thetaf from BC " << thetaf<< nl;

    // Visualization.
    auto& clangleBoundaryField = contactLineAngle_.boundaryFieldRef();
    auto& clanglePatchField = clangleBoundaryField[patchIndex];
    clanglePatchField *= 0.0;
    contactLineAngle_*=0.0;

    forAll(thetaf, faceI)
    {
       // Calculation of distance between cell centre and face centre
       const label cellID = faceOwner[patch.start() + faceI]; 
       vector cellCentre = C[cellID];
       const label faceID =  patch.start() + faceI; // Global face ID
       vector faceCentre = Cf[faceID];
       vector dw = cellCentre - faceCentre;
    
        if(hasContactLine(faceI))
        {
            //Info << nl << " From BC: thetaf " << thetaf[faceI] << " in cell " << cellID << nl;
            //Info << " and normal " << nHat[faceI] << nl;
            contactLineAngle_[faceOwner[patch.start() + faceI]] = thetaf[faceI]; 

            //Info << "Before: CellID : " << cellID << " ThetaD = " << thetaf[faceI] << " and uwall " << uwall[faceI]<< nl;
            // Visualize current PLIC contact angle as a cell-centered value.

            if
            (
               ( (thetaf[faceI] < thetaA_) && 
                (uwall[cellID] <0) )
                ||
                ( (thetaf[faceI] > thetaR_) && 
                (uwall[cellID] >0) )
            ) //advancing
            {
                //Info << " From BC : inside hysteresis" <<nl;
                    beta = 0;
                    //thetaf[faceI];
            }

            if((beta==1) && (uwall[faceI]!=0))
            {
                std::ofstream outFile("patchAlphaField.csv");
                // Check if the file is open
                if (outFile.is_open()) 
                {
                        outFile << 2024 << std::endl; // Writing each element on a new line
                }
                    // Close the file
                    outFile.close();
    
                //Info << " Inside beta==1" << nl;

                //thetaf[faceI] = theta0_;
                scalar uTheta_ = 1;
                thetaf[faceI] = theta0 + (thetaA_ - thetaR_)*tanh(uwall[faceI]/uTheta_);
                // thetaf[faceI] =  theta0_ * (1 - tanh(uwall[faceI]/uTheta_) * sign(uwall[faceI]/uTheta_) )
                //                 - thetaA_ * neg(uwall[faceI]) * tanh(uwall[faceI]/uTheta_)
                //                 + thetaR_ * pos(uwall[faceI]) * tanh(uwall[faceI]/uTheta_);
                // thetaf[faceI] =  min(
                //     180 / constant::mathematical::pi * (pow(ct_  * mag(Ca[faceI]) +
                //         pow(theta0_ * constant::mathematical::pi / 180, 3),
                //         0.3333333)),
                //     scalar(135)
                // );
            }
            //Info << nl << "After: CellID : " << cellID << " ThetaD = " << thetaf[faceI] << " and uwall " << uwall[faceI]<< nl;
            clanglePatchField[faceI] = thetaf[faceI];             
        }
    }
                
    // if(beta==0)
    // {
    //     // Info << " Inside beta=0" << nl;
    //     // forAll(pAlphaField, i)
    //     // {
    //     //     pAlphaField[i] = pAlphaFieldFile[i];
    //     // }
    //     // for (auto& vec : pUField) {
    //     //     for (double num : vec) {
    //     //         vec[num]*=scalar(0);
    //     //     }
    //     // }
        
    // }
    // alpha1_.write();  
    // U_.write();
    return thetafTmp;

}


void Foam::robinBC::write(Ostream& os) const
{
    alphaContactAngleTwoPhaseFvPatchScalarField::write(os);
    os.writeKeyword("theta0") << theta0_ << token::END_STATEMENT << nl;
    os.writeKeyword("ct") << ct_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        robinBC
    );
}


// ************************************************************************* //
