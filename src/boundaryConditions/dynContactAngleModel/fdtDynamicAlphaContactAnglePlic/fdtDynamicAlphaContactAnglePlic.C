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

\*---------------------------------------------------------------------------*/

#include "Identity.H"
#include "error.H"
#include "fdtDynamicAlphaContactAnglePlic.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "messageStream.H"
#include "scalarField.H"
#include "volFieldsFwd.H"
#include "volMesh.H"
#include "volFields.H"
#include "unitConversion.H"
#include "cutFacePLIC.H"
#include <tuple>
#include <vector>

// * * * * * * * * * * * Private Member Functions * * * * * * * * * * * * * * //

bool Foam::fdtDynamicAlphaContactAnglePlic::hasContactLine(label faceI) const
{
    // TODO: only works for 2D, remove
    //const auto& connectedFaceLabels = this->patch().patch().faceFaces()[faceI];
    //const auto alphaInternalTmp = this->patchInternalField();
    //const auto& alphaInternal = alphaInternalTmp.cref();
    //const scalar tol = 1.0e-8;

    //forAll(connectedFaceLabels, I)
    //{
        //// Don't use zero, but a small value as tolerance
        //if ((alphaInternal[connectedFaceLabels[I]] < tol) && (alphaInternal[faceI] > tol))
        //{n einer Kapillare, oder bei einer Ausbreitung eines Tropfens
            //contactLine = true;
            //Pout<< "local face ID " << faceI
                //<< "; connected face ID " << connectedFaceLabels[I]
                //<< endl;
        //}
    //}

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

std::tuple<Foam::point,Foam::point> Foam::fdtDynamicAlphaContactAnglePlic::closestPoints
(
    point pA, vector dirA, point pB, vector dirB
) const 
{
    // Ensure lines are not parallel
    if ((1.0 - mag((dirA & dirB))/(mag(dirA)*mag(dirB))) < SMALL)
    {
        // Lines are parallel, throw error or something
    }

    // The connection between the closest points on two different lines
    // is a line that is perpendicular to both lines. This fact is used to
    // construct a 2x2 system A*x=b with scalar factors to the direction
    // vector as unknowns.
    
    // Right hand side
    scalar b1 = (pA - pB) & dirA;
    scalar b2 = (pA - pB) & dirB;
    scalar a11 =  dirA & dirA;
    scalar a12 = -dirB & dirA;
    scalar a21 =  dirA & dirB;
    scalar a22 = -dirB & dirB;
    scalar detA = a11*a22 - a12*a21;
    
    // Use analytic formula for inverse to compute the sought after
    // scalar factors
    scalar lambdaA = (a22*b1 - a12*b2)/detA;
    scalar lambdaB = (-a21*b1 + a11*b2)/detA;

    point closestOnA = pA + lambdaA*dirA;
    point closestOnB = pB + lambdaB*dirB;

    return std::tuple<point,point>{closestOnA,closestOnB};
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fdtDynamicAlphaContactAnglePlic::
fdtDynamicAlphaContactAnglePlic
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    alphaContactAngleTwoPhaseFvPatchScalarField(p, iF),
    theta0_(0.0),
    thetaA_(0.0),
    thetaR_(0.0),
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


Foam::fdtDynamicAlphaContactAnglePlic::
fdtDynamicAlphaContactAnglePlic
(
    const fdtDynamicAlphaContactAnglePlic& gcpsf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    alphaContactAngleTwoPhaseFvPatchScalarField(gcpsf, p, iF, mapper),
    theta0_(gcpsf.theta0_),
    thetaA_(gcpsf.thetaA_),
    thetaR_(gcpsf.thetaR_),
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


Foam::fdtDynamicAlphaContactAnglePlic::
fdtDynamicAlphaContactAnglePlic
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    alphaContactAngleTwoPhaseFvPatchScalarField(p, iF, dict),
    theta0_(dict.get<scalar>("theta0")),
    thetaA_(dict.get<scalar>("thetaA")),
    thetaR_(dict.get<scalar>("thetaR")),
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


Foam::fdtDynamicAlphaContactAnglePlic::
fdtDynamicAlphaContactAnglePlic
(
    const fdtDynamicAlphaContactAnglePlic& gcpsf
)
:
    alphaContactAngleTwoPhaseFvPatchScalarField(gcpsf),
    theta0_(gcpsf.theta0_),
    thetaA_(gcpsf.thetaA_),
    thetaR_(gcpsf.thetaR_),
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


Foam::fdtDynamicAlphaContactAnglePlic::
fdtDynamicAlphaContactAnglePlic
(
    const fdtDynamicAlphaContactAnglePlic& gcpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    alphaContactAngleTwoPhaseFvPatchScalarField(gcpsf, iF),
    theta0_(gcpsf.theta0_),
    thetaA_(gcpsf.thetaA_),
    thetaR_(gcpsf.thetaR_),
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
std::vector<std::tuple<Foam::point,Foam::vector>>
Foam::fdtDynamicAlphaContactAnglePlic::interpolatePLICNormals
(
    label faceI,
    const volVectorField& interfaceCentre,
    const volVectorField& interfaceNormal
) const
{
    std::vector<std::tuple<point,vector>> pointsAndNormals{};

    const auto& mesh = interfaceCentre.mesh();
    auto cellI = mesh.faceOwner()[faceI];

    // TT: First, only try the face neighbour cells and see if this works
    // reliably. If not, we may need also the edge and point neighbouring cells
    const auto& cellToCell = mesh.cellCells();
    const auto& cells = mesh.cells();
    const auto& thisCell = cells[cellI];

    for(auto cellK : cellToCell[cellI])
    {
        // TODO (TT): when initializing normals from the ones given by SMCIA,
        // we have non-zero normals in every cell. We need an additional check to
        // rule them out as candidates.
        if (magSqr(interfaceNormal[cellK]) != 0.0)
        {
            // Get the face connecting cellI and cellK
            label connectingFaceI = -1;
            const auto& otherCell = cells[cellK];
            for (auto thisFaceI : thisCell)
            {
                for (auto otherFaceK : otherCell)
                {
                    if (thisFaceI == otherFaceK)
                    {
                        connectingFaceI = thisFaceI;
                        break;
                    }
                }

                if (connectingFaceI != -1)
                {
                    break;
                }
            }
            
            // Get the cuts of the face connecting cellI and cellK with their
            // respective PLICs.
            // Only continue with interpolation in case both PLICs intersect the
            // face shared by the two cells

            // Note (TT): the field 'interfaceNormals' does not provide unit normals,
            // but probably normals scaled with the area of the corresponding PLIC
            auto normalI = interfaceNormal[cellI];
            normalI.normalise();
            auto normalK = interfaceNormal[cellK];
            normalK.normalise();

            cutFacePLIC cutFaceI(mesh);
            label cutStatusI = cutFaceI.calcSubFace
                (
                    connectingFaceI,
                    normalI,
                    interfaceCentre[cellI]
                );

            // Debug
            //Info<< "Cut status of cellI: " << cutStatusI << endl;
            
            cutFacePLIC cutFaceK(mesh);
            label cutStatusK = cutFaceK.calcSubFace
                (
                    connectingFaceI,
                    normalK,
                    interfaceCentre[cellK]
                );
            // Debug
            //Info<< "Cut status of cellK: " << cutStatusK << endl;
            //Info<< "Normal_I: " << interfaceNormal[cellI] << ", centroid_I: " << interfaceCentre[cellI]
            //    << ", PLIC area_I: " << mag(interfaceNormal[cellI]) << nl
            //    << "Normal_K: " << interfaceNormal[cellK] << ", centroid_K: " << interfaceCentre[cellK]
            //    << ", PLIC area_K: " << mag(interfaceNormal[cellK]) << endl;
            if (cutStatusK != 0 || cutStatusI != 0)
            {
                continue;
            }

            // Perform interpolation. The idea to find the proper point to interpolate to
            // is sketeched below:
            //  1) The PLIC cuts of the face do not coincide in general. Thus, compute an
            //      average intersection line L.
            //  2) Compute the line connecting both PLIC centres CC.
            //  3) On L, find the point closest to CC. This is the point PI to interpolate to.
            //  4) As interpolation weights the distances between the PLIC centres and PI
            //     are used.
            const auto& cutPointsI = cutFaceI.surfacePoints();
            const auto& cutPointsL = cutFaceK.surfacePoints();

            // Debug
            //Info<< "Cut points of face I: " << cutPointsI << endl;
            //Info<< "Cut points of face K: " << cutPointsL << endl;

            point LA{0,0,0};
            point LB{0,0,0};

            // Average the points that are closer to each other
            if (magSqr(cutPointsI[0] - cutPointsL[0]) < magSqr(cutPointsI[0] - cutPointsL[1]))
            {
                LA = 0.5*(cutPointsI[0] + cutPointsL[0]);
                LB = 0.5*(cutPointsI[1] + cutPointsL[1]);
            }
            else
            {
                LA = 0.5*(cutPointsI[0] + cutPointsL[1]);
                LB = 0.5*(cutPointsI[1] + cutPointsL[0]);
            }

            // Construct direction of line L
            // Debug
            //Info<< "Averaged points from cutting mesh face with PLICs: "
            //    << "LA = " << LA << ", LB = " << LB << endl; 
            vector ldir = LB - LA;

            // Construct line between PLIC centroids
            point CA{interfaceCentre[cellI]};
            vector cdir = interfaceCentre[cellK] - interfaceCentre[cellI];
            // Debug
            //Info<< "Vector connecting PLIC centroids: " << cdir << endl;

            //auto [pInterp, noUsed] = closestPoints(LA, ldir, CA, cdir);
            point pInterp{}, notUsed{};
            std::tie(pInterp, notUsed) = closestPoints(LA, ldir, CA, cdir);

            // Debug
            //Info<< "Interpolating PLIC normals to " << pInterp << endl;

            // Use distances to pInterp as interpolation weights
            auto distI = mag(interfaceCentre[cellI] - pInterp);
            auto distK = mag(interfaceCentre[cellK] - pInterp);

            vector nInterp = (distI*normalK + distK*normalI)/(distI + distK);
            nInterp.normalise();

            // Debug
            // Plausibility check: the angle enclosed between the the two PLIC normals should
            //  be larger than the angle enclosed by either of them with the interpolated normal.
            //  For simplicity, it is assumed that all enclosed angles are smaller than 90 degrees.
            //auto plicAngle = radToDeg(Foam::acos(normalI & normalK));
            //auto angleIInterp = radToDeg(Foam::acos(normalI & nInterp));
            //auto angleKInterp = radToDeg(Foam::acos(normalK & nInterp));
            //Info<< "Angle enclosed by PLIC normals: " << plicAngle << nl
            //    << "Angle between PLIC I and interpolated normal: " << angleIInterp << nl
            //    << "Angle between PLIC K and interpolated normal: " << angleKInterp << nl
            //    << "Sum of interpolated angles: " << angleIInterp + angleKInterp
            //    << endl;

            pointsAndNormals.push_back(std::make_tuple(pInterp, nInterp));
        }
    }

    return pointsAndNormals;
}

Foam::vector Foam::fdtDynamicAlphaContactAnglePlic::extrapolatePLICNormals
(
    point targetPoint,
    point plicCentre,
    vector plicNormal,
    const std::vector<std::tuple<point,vector>>& pointsAndNormals
) const
{
    // Ensure we have a unit normal
    plicNormal.normalise();

    // Construct the plane spanned by the plic centroid, its normal and the
    // line to the target point
    auto planeNormal = (targetPoint - plicCentre) ^ plicNormal;
    planeNormal.normalise();

    // Project the points and normals into the plane
    std::vector<point> projectedPoints(pointsAndNormals.size());
    std::vector<vector> projectedNormals(pointsAndNormals.size());
    auto projector = (tensor::I - (planeNormal*planeNormal));

    //Info<< "Projecting points and normals" << endl;
    for (uLabel I = 0; I != pointsAndNormals.size(); ++I)
    {
        point p{};
        vector n{};

        std::tie(p, n) = pointsAndNormals[I];

        projectedPoints[I] = projector & (p - plicCentre);
        projectedNormals[I] = projector & n;
        projectedNormals[I].normalise();
    }

    // Transform points and normals into the local coordinate system
    // The local, 2D coordinate system is defined as follows:
    //  - The connection from the PLIC centroid to the target point
    //      defines the positive x-direction.
    //  - The PLIC normal defines the positive y-direction.
    std::vector<point> transformedPoints(pointsAndNormals.size());
    std::vector<vector> transformedNormals(pointsAndNormals.size());

    // The transformation into the local coordinate system can be done
    // by scalar products with the vectors representing base vectors
    // of the local coordinate system.
    vector exl = (targetPoint - plicCentre);
    exl.normalise();
    vector eyl = plicNormal;
    // The z component in the local coordinate system is zero by contruction
    // due to the projection into a plane done above.

    for (uLabel I = 0; I != pointsAndNormals.size(); ++I)
    {
        point p = projectedPoints[I];
        transformedPoints[I] = point{p&exl, p&eyl, 0.0};

        vector n = projectedNormals[I];
        transformedNormals[I] = vector{n&exl, n&eyl, 0.0};
    }

    // Construct and solve the least-squares fit for a parabola
    // TODO: document this somewhere, noone can follow this without derivation
    // The linear least squares problem using the normal formulation results
    // in a 2x2 system to fit the parabola. The third coefficient of the parabola
    // is set to zero per construction so the parabola passes through the PLIC centroid.
    // linCoeff is the first unknown of the system, quadCoeff the second one,
    // according to 
    // f(x) = linCoeff*x + quadCoeff*x^2

    // The contributions from the PLIC centroid and its normals are not included in
    // pointsAndNormals, so they are used here to initialize the values.
    scalar A11 = 1.0;
    scalar A12 = 0.0;
    scalar A22 = 0.0;
    scalar b1 = 0.0;
    scalar b2 = 0.0;

    for (uLabel I = 0; I != pointsAndNormals.size(); ++I)
    {
        auto p = transformedPoints[I];
        auto n = transformedNormals[I];
        A11 += pow(p.x(), 2) + pow(n.y(), 2);
        A12 += pow(p.x(), 3) + 2*p.x()*pow(n.y(), 2);
        A22 += pow(p.x(), 4) + 4*pow(p.x(), 2)*pow(n.y(), 2);

        b1 += p.x()*p.y() - p.x()*n.x()*n.y();
        b2 += pow(p.x(), 2)*p.y() - 2*pow(p.x(), 2)*n.x()*n.y();
    }

    scalar A21 = A12;
    scalar detA = A11*A22 - A21*A12;

    scalar linCoeff = 1.0/detA*(A22*b1 - A21*b2);
    scalar quadCoeff = 1.0/detA*(-A12*b1 + A11*b2);

    // Per construction of the local coordinate system, the x-value of the target
    // point is the distance between it and the PLIC centroid
    auto xtp = mag(targetPoint - plicCentre);

    // Compute the interface normal at the target point, transform it
    // to the global coordinate system and return it.
    // Transformation from the local to the global coordinate system is done
    // by computing the product of the vector component in local coordinate
    // system with the corresponding base vector expressed in the global
    // coordinate system.

    // x-component of the normal in local coordinate system is the derivative of
    // the parabola at xtp
    auto targetNormal_X_local = linCoeff + 2*quadCoeff*xtp;
    auto targetNormal_Y_local = -xtp;

    auto targetNormal = targetNormal_X_local*exl + targetNormal_Y_local*plicNormal; 
    targetNormal.normalise();

    // Check orientation of normal with respect to PLIC normal, they should enclose
    // an angle smaller than 90 degree
    if ((targetNormal & plicNormal) < 0.0)
    {
        targetNormal *= -1.0;
    }

    // Debug
    //Info<< "PLIC normal: " << plicNormal << ", extraploted normal: "
    //    << targetNormal << endl;

    return targetNormal;
}

Foam::tmp<Foam::scalarField>
Foam::fdtDynamicAlphaContactAnglePlic::theta
(
    const fvPatchVectorField& Up,
    const fvsPatchVectorField& nHat 
) const
{
    // Get the interface normals at the wall.
    const vectorField nf(patch().nf());

    // Find the direction of the interface parallel to the wall
    vectorField nWall(nHat - (nf & nHat)*nf);
    // Normalise nWall
    nWall /= (mag(nWall) + SMALL);

    // Calculate contact line velocity relative to the wall 
    // velocity for moving walls. 
    vectorField Uwall(Up.patchInternalField() - Up);

    // Calculate component of the contact line velocity Uwal in the direction
    // of the interface normal tagential to the wall.
    scalarField uwall(nWall & Uwall);

    // Fetch physical properties
    // TODO(TM): calculate single-field properties?
    // TODO(TM): check if \rho1\nu1 > \rho2\nu2

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
    scalarField CaWall(muwall*uwall/sigmap.value());

    // Compute the contact angles at the wall.
    tmp<scalarField> thetafTmp = Foam::radToDeg(Foam::acos(nHat & nf));
    scalarField& thetaf = thetafTmp.ref();

    // Visualization.
    const auto& patch = this->patch();
    const label patchIndex =patch.index();
    auto& clangleBoundaryField = contactLineAngle_.boundaryFieldRef();
    auto& clanglePatchField = clangleBoundaryField[patchIndex];

    // Reset clangePatchField to reset contact angles from previous time steps
    clanglePatchField = 0.0;
    
    // Visualization of the contact angle
    const fvMesh& mesh = nu.mesh(); 

    // Contact line length calculation data
    const auto& magSf = mesh.magSf();
    const auto& boundaryMagSf = magSf.boundaryField();
    const auto& patchMagSf = boundaryMagSf[patchIndex];

    const auto& dCoeffs = mesh.deltaCoeffs();  
    const auto& boundaryDcoeffs = dCoeffs.boundaryField(); 
    const auto& patchDcoeffs = boundaryDcoeffs[patchIndex];

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
    const auto alphaName = IOobject::groupName
    (
        "alpha", 
        this->internalField().group()
    );

    bool hasNormals = db.found(normalsName);
    if (!hasNormals)
    {
        // This BC is updated before interface reconstruction.
        // Do nothing if PLIC fields are not available in the registry. 
        thetaf = GREAT;
        return thetafTmp;
    }

    bool hasCentres = db.found(centresName);
    if (!hasCentres)
    {
        // This BC is updated before interface reconstruction.
        // Do nothing if PLIC fields are not available in the registry. 
        thetaf = GREAT; 
        return thetafTmp;
    }

    const volVectorField& interfaceNormal = 
        db.lookupObject<volVectorField>(normalsName);

    const volVectorField& interfaceCentre = 
        db.lookupObject<volVectorField>(centresName);

    const volScalarField& alpha = 
        db.lookupObject<volScalarField>(alphaName);

    cutFacePLIC cutFace(mesh);

    // For all boundary faces
    //Info << "\n--- Entering loop over FDT boundary faces ---" << endl;

    // Summarized debug info
    uLabel nCLcells = 0;
    scalar caMin = 180.0;
    scalar caMax = 0.0;
    scalar caMean = 0.0;
    uLabel nHysteresisCells = 0;
    scalar dThetaMin = 180.0;
    scalar dThetaMax = -180.0;
    scalar dThetaMean = 0.0;

    forAll(thetaf, faceI)
    {
        const label cellI = patch.faceCells()[faceI];

        // If we are in a contact-line cell
        // TODO: at least in the first time step, checking for a non-zero normal is not sufficient.
        //      The initialization sets the interface normals in every cell.
        //      Thus, we need to check the alpha field or something else in addition.
        if (mag(nHat[faceI]) > 0 && (alpha[cellI] > 1.0e-8) && (alpha[cellI] < (1.0-1.0e-8))) //TODO(TM): && hasContactLine(faceI))
        {
            const label globalFaceI = patch.start() + faceI;

            // Extrapolate PLIC normals to contact line by first interpolating to the PLIC from
            // neighbours and then do the actual extrapolation
            if (hasContactLine(faceI))
            {
                auto pointsAndNormals = interpolatePLICNormals(globalFaceI, interfaceCentre, interfaceNormal);
                //Info<< "Found " << pointsAndNormals.size() << " interpolated normals." << endl;

                // Debug: catch the case that against expectation no neighbouring cells with
                // PLIC normals are found
                if (pointsAndNormals.size() > 0)
                {
                    cutFacePLIC boundaryCut(mesh);
                    auto normalI = interfaceNormal[cellI];
                    normalI.normalise();
                    boundaryCut.calcSubFace(globalFaceI, normalI, interfaceCentre[cellI]);
                    auto intersections = boundaryCut.surfacePoints();
                    auto contactLineCentre = 0.5*(intersections[0] + intersections[1]);

                    //Info<< "Extrapolating normal to contact line" << endl;

                    auto contactLineNormal = extrapolatePLICNormals
                                                   (
                                                        contactLineCentre,
                                                        interfaceCentre[cellI],
                                                        normalI,
                                                        pointsAndNormals
                                                    );
                    auto contactAngle = radToDeg(Foam::acos(contactLineNormal & nf[faceI]));
                    //Info<< "Contact angle difference between PLIC and extrapolated normal: " << mag(thetaf[faceI] - contactAngle) << endl;
                    thetaf[faceI] = contactAngle;
                }
                else
                {
                    // Throw error
                    //Info<< "Warning: did not find ant normals to interpolate for cell " << cellI << endl;
                }
            }

            // Visualize current PLIC contact angle as a cell-centered value.
            contactLineAngle_[cellI] = thetaf[faceI];

            //Pout << "theta_old = " << thetaf[faceI] << endl;
            if (thetaf[faceI] < thetaR_) // Receding regime
            {
                thetaf[faceI] = thetaR_;
            }
            else if (thetaf[faceI] > thetaA_) // Advancing regime
            {
                thetaf[faceI] = thetaA_;
            }
            else // Hysteresis regime
            {
                // Equation 32 in the manuscript.
                scalar Cstar = patchMagSf[faceI] * patchDcoeffs[faceI] *
                    (thetaA_ - thetaR_) * muwall[faceI]  / 
                    (
                         mag(Foam::cos(Foam::degToRad(thetaA_)) - 
                             Foam::cos(Foam::degToRad(thetaR_))) * sigmap.value()
                    );

                // Cut the face.
                label cutStatus = cutFace.calcSubFace
                (
                    faceI,
                    interfaceNormal[cellI],
                    interfaceCentre[cellI]
                );

                // If the face is actually cut
                if(cutStatus == 0)
                {
                    const auto& cutPoints = cutFace.surfacePoints(); 
                    // Geometrical contact line length
                    scalar Cll = Foam::mag (cutPoints[1] - cutPoints[0]);  
                    // Scale Cstar with Cll 
                    Cstar = Cstar / max(Cll, 0.1*Foam::sqrt(patchMagSf[faceI]));

                } else 
                {
                    // Use largest possible contact line length sqrt(magSf_b)
                    Cstar = Cstar / Foam::sqrt(patchMagSf[faceI]);
                }
                // Equation 31 in the manuscript. Limit change to 5.0 degrees maximum
                scalar dtheta = min(5.0, Cstar * mag(uwall[faceI]));
                // scalar dtheta = Cstar * mag(uwall[faceI]);

                if (uwall[faceI] < 0)
                {
                    thetaf[faceI] += dtheta;
                    //Pout << "Hysteresis mode advancing, " 
                    //    << " dtheta = " << dtheta 
                    //    << " uwall = " << uwall[faceI] 
                    //    << endl;
                }
                else if (uwall[faceI] > 0)
                {
                    thetaf[faceI] -= dtheta;
                    //Pout << "Hysteresis mode receding, " 
                    //    << " dtheta = " << -dtheta 
                    //    << " uwall = " << uwall[faceI] 
                    //    << endl;
                }
                else
                {
                    //Pout << "Do nothing thetaf = " << thetaf[faceI] << endl;
                }

                // Collect summarized debug info
                ++nHysteresisCells;
                dThetaMin = min(dThetaMin, dtheta);
                dThetaMax = max(dThetaMax, dtheta);
                dThetaMean += dtheta;
            }
            // Visualize the new dynamic contact angle as a face-centered value.
            clanglePatchField[faceI] = thetaf[faceI];
            //Pout << "Contact line on face " << faceI
            //     << "Cell ID " << nu.mesh().faceOwner()[faceI + this->patch().start()]
            //     << "\n\ttheta = " << thetaf[faceI]
            //     << "\n\tnWall = " << nWall[faceI]
            //     << "\n\tuwall = " << uwall[faceI]
            //     <<endl;
            
            // Collect info for debug summary
            ++nCLcells;
            caMean += thetaf[faceI];
            caMin = min(thetaf[faceI], caMin);
            caMax = max(thetaf[faceI], caMax);
        }
    }

    caMean /= nCLcells;

    // Output debug summary
    Pout<< "Contact line info:" << nl
        << "\tnumber of contact line cells: " << nCLcells << nl
        << "\tcontact angles (min/mean/max): "
        << caMin << ",\t" << caMean << ",\t" << caMax
        << endl;
    if (nHysteresisCells > 0)
    {
        Pout<< "\tnumber of hysteresis cells: " << nHysteresisCells << nl
            << "\tdtheta values (min/mean/max): "
            << dThetaMin << ",\t" << dThetaMean/nHysteresisCells << ",\t"
            << dThetaMax << endl;
    }

    return thetafTmp; 
}


void Foam::fdtDynamicAlphaContactAnglePlic::write(Ostream& os) const
{
    alphaContactAngleTwoPhaseFvPatchScalarField::write(os);
    os.writeEntry("theta0", theta0_);
    os.writeEntry("thetaA", thetaA_);
    os.writeEntry("thetaR", thetaR_);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        fdtDynamicAlphaContactAnglePlic
    );
}


// ************************************************************************* //
