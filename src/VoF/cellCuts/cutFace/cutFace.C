/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 DHI
    Copyright (C) 2018-2019 Johan Roenby
    Copyright (C) 2019-2020 DLR
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

#include "cutFace.H"
#include <vector>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

int Foam::cutFace::debug = 0;
std::vector<double> Foam::cutFace::pOld = {0.0, 0.0, 0.0};
std::unordered_map<Foam::label, std::vector<Foam::point>> Foam::cutFace::faceIntersectionsMap;
int Foam::cutFace::timeIndex_ = -100;
// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * *

void Foam::cutFace::storeAndCalculateIntersections(
    label faceI,
    const scalarList& pointStatus,
    const pointField& points,
    const face& f,
    DynamicList<point>& subFacePoints,
    DynamicList<point>& surfacePoints,
    label firstFullySubmergedPoint
)
{
    std::vector<point> intersectionPoints;

    // Loop starting at firstFullySubmergedPoint
    for (label i = firstFullySubmergedPoint; i < firstFullySubmergedPoint + f.size(); ++i)
    {
        label idx = i % f.size();
        label nextIdx = (i + 1) % f.size();

        // Append fluid points
        if (pointStatus[idx] > 0)
        {
            subFacePoints.append(points[f[idx]]);
        }
        // Append cut points
        else if (pointStatus[idx] == 0)
        {
            subFacePoints.append(points[f[idx]]);
            surfacePoints.append(points[f[idx]]);
        }

        // Calculate intersection points for edges
        if (
            (pointStatus[idx] < 0 && pointStatus[nextIdx] > 0) ||
            (pointStatus[idx] > 0 && pointStatus[nextIdx] < 0))
        {
            // Compute intersection point on the edge
            label nextP = f.nextLabel(idx);
            vector dir = points[nextP] - points[f[idx]];
            scalar weight = (0.0 - pointStatus[idx]) /
                            (pointStatus[nextIdx] - pointStatus[idx]);

            point intersection = points[f[idx]] + weight * dir;

            // Store the intersection point
            intersectionPoints.push_back(intersection);

            // Append intersection points to subFacePoints and surfacePoints
            subFacePoints.append(intersection);
            surfacePoints.append(intersection);
        }
    }

    // Store the intersections in the static map
    faceIntersectionsMap[faceI] = intersectionPoints;
}


void Foam::cutFace::calcSubFace
(
    const label faceI,
    const scalarList& pointStatus,
    label firstFullySubmergedPoint,
    DynamicList<point>& subFacePoints,
    DynamicList<point>& surfacePoints,
    label& faceStatus,
    vector& subFaceCentre,
    vector& subFaceArea
)
{
    const pointField& points = mesh_.points();
    const face& f = mesh_.faces()[faceI];

    if (firstFullySubmergedPoint == -1) // In gas phase
    {
        faceStatus = 1;
        subFaceCentre = Zero;
        subFaceArea = Zero;
        return;
    }

    // Call the function to calculate and store intersections
    storeAndCalculateIntersections(faceI, pointStatus, points, f, subFacePoints, surfacePoints, firstFullySubmergedPoint);

    // Retrieve stored intersections (optional, for debugging or post-processing)
    const std::vector<point>& intersections = faceIntersectionsMap[faceI];

    // Additional logic for sub-face area and center
    if (subFacePoints.size() >= 3)
    {
        faceStatus = 0;
        calcSubFaceCentreAndArea(subFacePoints, subFaceCentre, subFaceArea);
    }
    else
    {
        faceStatus = -1;
    }







    // label patchi = 0;

    // for (label patchID = 0; patchID < mesh_.boundaryMesh().size(); ++patchID)
    // {
    //     const auto patch = mesh_.boundaryMesh()[patchID];
    //     const label patchStart = patch.start();
    //     const label patchEnd = patchStart + patch.size();

    //     // Check if the face ID falls within this patch's range
    //     if (faceI >= patchStart && faceI < patchEnd)
    //     {
    //         patchi = patchID;
    //     }
    // }

    // const polyBoundaryMesh& boundaryMesh = mesh_.boundaryMesh();
    // const label start = boundaryMesh[patchi].start();

    // // if (!mesh_.isInternalFace(faceI) && mesh_.boundaryMesh()[patchi].name() =="bottom")
    // //     Info << "\nStart cutFace-1 -> calcSubFace:\t Patch:\t" << mesh_.boundaryMesh()[patchi].name() << "\nhas the global face:\t" << faceI << "\nand local face:\t" << faceI-start << "\n";
    // const pointField& points = mesh_.points();
    // const face& f = mesh_.faces()[faceI];

    // if (firstFullySubmergedPoint == -1) // is in gasPhase
    // {
    //     faceStatus = 1;
    //     subFaceCentre = Zero;
    //     subFaceArea = Zero;
    //     return;
    // }

    // // loop face and append the cuts
    // // loop starts at firstFullySubmergedPoint
    // for
    // (
    //     label i = firstFullySubmergedPoint;
    //     i < firstFullySubmergedPoint + f.size();
    //     ++i
    // )
    // {
    //     // max two points are appended during one cycle
    //     label idx = i % f.size();
    //     label nextIdx = (i + 1) % f.size();

    //     if (pointStatus[idx] > 0) // append fluid point
    //     {
    //         subFacePoints.append(points[f[idx]]);
    //     }
    //     else if (pointStatus[idx] == 0) // append cut point
    //     {
    //         subFacePoints.append(points[f[idx]]);
    //         surfacePoints.append(points[f[idx]]);
    //     }

    //     if
    //     (
    //         (pointStatus[idx] < 0 && pointStatus[nextIdx] > 0) ||
    //         (pointStatus[idx] > 0 && pointStatus[nextIdx] < 0)
    //     ) // cut on edge cut Value is zero
    //     {
    //         label nextP = f.nextLabel(idx);
    //         vector dir = points[nextP] - points[f[idx]];
    //         scalar weight =
    //             (0.0 - pointStatus[idx]) /
    //             (pointStatus[nextIdx] - pointStatus[idx]); // cutValue is zero
    //         if (!mesh_.isInternalFace(faceI) && mesh_.boundaryMesh()[patchi].name() =="bottom" && faceI==208 && pOld[1]==0.0)
    //         {
    //             Info << "\nFace " << faceI << " Debugging for p:\n";
    //             Info << "Before:\t " << pOld << endl;
    //         }
    //         point p = points[f[idx]] + weight * dir; // intersection point
    //         if (!mesh_.isInternalFace(faceI) && mesh_.boundaryMesh()[patchi].name() =="bottom" && faceI==208 && p.y()==0.0)
    //         {
    //             Info << "after:\t " << p << " and pold " << pOld << endl;
    //         }
    //         Info << "timeIndex " << timeIndex_ << " mesh time " <<mesh_.time().timeIndex() << "\n";
    //         if (!mesh_.isInternalFace(faceI) && mesh_.boundaryMesh()[patchi].name() =="bottom" && faceI==208 && p.y()==0.0)
    //         {
    //             if(timeIndex_==-100)
    //             {
    //                 Info << "In first loop\n";
    //                 pOld = {p.x(), p.y(), p.z()};
    //                 timeIndex_= mesh_.time().timeIndex();
    //             }
    //             if(timeIndex_ <= mesh_.time().timeIndex())
    //             {
    //                 Info << "In second loop\n";
    //                 timeIndex_= mesh_.time().timeIndex();
    //                 p = point(pOld[0], pOld[1], pOld[2]);
    //             }
                
    //             Info << "Intersection point of interface and face " << faceI << " is " << p << endl;
    //         }
    //         if (!mesh_.isInternalFace(faceI) && mesh_.boundaryMesh()[patchi].name() =="bottom" && faceI==208 && p.y()==0.0)
    //         {
    //             Info << "aftersquared:\t " << p << " and pold " << pOld<< endl;
    //         }
    //         subFacePoints.append(p);
    //         surfacePoints.append(p);
    //     }
    // }

    // if (subFacePoints.size() >= 3)
    // {
    //     faceStatus = 0;
    //     calcSubFaceCentreAndArea(subFacePoints, subFaceCentre, subFaceArea);
    // }
    // else
    // {
    //     faceStatus = -1;
    // }


    // if (!mesh_.isInternalFace(faceI) && mesh_.boundaryMesh()[patchi].name() =="bottom")
    //     Info << "\nEnd cutFace-1 -> calcSubFace:\t Patch:\t" << mesh_.boundaryMesh()[patchi].name() << "\nhas the global face:\t" << faceI << "\nand local face:\t" << faceI-start << "\n";
}


void Foam::cutFace::calcSubFace
(
    const label faceI,
    const scalarList& pointStatus,
    const scalarList& weights,
    label firstFullySubmergedPoint,
    DynamicList<point>& subFacePoints,
    DynamicList<point>& surfacePoints,
    label& faceStatus,
    vector& subFaceCentre,
    vector& subFaceArea
)
{
    label patchi = 0;

    for (label patchID = 0; patchID < mesh_.boundaryMesh().size(); ++patchID)
    {
        const auto patch = mesh_.boundaryMesh()[patchID];
        const label patchStart = patch.start();
        const label patchEnd = patchStart + patch.size();

        // Check if the face ID falls within this patch's range
        if (faceI >= patchStart && faceI < patchEnd)
        {
            patchi = patchID;
        }
    }

    const polyBoundaryMesh& boundaryMesh = mesh_.boundaryMesh();
    const label start = boundaryMesh[patchi].start();

    // if (!mesh_.isInternalFace(faceI) && mesh_.boundaryMesh()[patchi].name() =="bottom")
    //     Info << "\nStart cutFace-2 -> calcSubFace:\t Patch:\t" << mesh_.boundaryMesh()[patchi].name() << "\nhas the global face:\t" << faceI << "\nand local face:\t" << faceI-start << "\n";
    const pointField& points = mesh_.points();
    const face& f = mesh_.faces()[faceI];

    if (firstFullySubmergedPoint == -1) // is in gasPhase
    {
        faceStatus = 1;
        subFaceCentre = Zero;
        subFaceArea = Zero;
        return;
    }

    // loop face and append the cuts
    // loop starts at firstFullySubmergedPoint
    for
    (
        label i = firstFullySubmergedPoint;
        i < firstFullySubmergedPoint + f.size();
        ++i
    )
    {
        // max two points are appended during one cycle
        label idx = i % f.size();
        label nextIdx = (i + 1) % f.size();

        if (pointStatus[idx] > 0) // append fluid point
        {
            subFacePoints.append(points[f[idx]]);
        }
        else if (pointStatus[idx] == 0) // append cut point
        {
            subFacePoints.append(points[f[idx]]);
            surfacePoints.append(points[f[idx]]);
        }

        if
        (
            (pointStatus[idx] < 0 && pointStatus[nextIdx] > 0) ||
            (pointStatus[idx] > 0 && pointStatus[nextIdx] < 0)
        ) // cut on edge cut Value is zero
        {
            label nextP = f.nextLabel(idx);
            vector dir = points[nextP] - points[f[idx]];

            point p = points[f[idx]] + weights[idx] * dir;

            subFacePoints.append(p);
            surfacePoints.append(p);
        }
    }

    if (subFacePoints.size() >= 3)
    {
        faceStatus = 0;
        calcSubFaceCentreAndArea(subFacePoints, subFaceCentre, subFaceArea);
    }
    else
    {
        faceStatus = -1;
    }


    // if (!mesh_.isInternalFace(faceI) && mesh_.boundaryMesh()[patchi].name() =="bottom")
    //     Info << "\nEnd cutFace-2 -> calcSubFace:\t Patch:\t" << mesh_.boundaryMesh()[patchi].name() << "\nhas the global face:\t" << faceI << "\nand local face:\t" << faceI-start << "\n";
}


void Foam::cutFace::calcSubFace
(
    const face& f,
    const pointField& points,
    const scalarList& pointStatus,
    label firstFullySubmergedPoint,
    DynamicList<point>& subFacePoints,
    DynamicList<point>& surfacePoints,
    label& faceStatus,
    vector& subFaceCentre,
    vector& subFaceArea
)
{
    
    Info << "I am in cutFace -3 start\n";
    if (firstFullySubmergedPoint == -1) // in Gas
    {
        faceStatus = 1;
        subFaceCentre = Zero;
        subFaceArea = Zero;
        return;
    }

    // loop face and append the cuts
    for
    (
        label i = firstFullySubmergedPoint;
        i < firstFullySubmergedPoint + f.size();
        ++i
    )
    {
        // max two points are appended during one cycle
        label idx = i % f.size();
        label nextIdx = (i + 1) % f.size();

        if (pointStatus[idx] > 0) // append fluid point
        {
            subFacePoints.append(points[f[idx]]);
        }
        else if (pointStatus[idx] == 0) // append cut point
        {
            subFacePoints.append(points[f[idx]]);
            surfacePoints.append(points[f[idx]]);
        }

        if
        (
            (pointStatus[idx] < 0 && pointStatus[nextIdx] > 0) ||
            (pointStatus[idx] > 0 &&  pointStatus[nextIdx] < 0)
        )
        {
            label nextP = f.nextLabel(idx);
            vector dir = points[nextP] - points[f[idx]];
            scalar weight =
                (0.0 - pointStatus[idx]) /
                (pointStatus[nextIdx] - pointStatus[idx]);

            point p = points[f[idx]] + weight * dir;

            subFacePoints.append(p);
            surfacePoints.append(p);
        }
    }

    if (subFacePoints.size() >= 3)
    {
        faceStatus = 0;
        calcSubFaceCentreAndArea(subFacePoints, subFaceCentre, subFaceArea);
    }
    else
    {
        faceStatus = -1;
    }
        Info << "I am in cutFace -3 end \n";
}


void Foam::cutFace::calcSubFaceCentreAndArea
(
    DynamicList<point>& subFacePoints,
    vector& subFaceCentre,
    vector& subFaceArea
)
{
    const label nPoints = subFacePoints.size();

    // If the face is a triangle, do a direct calculation for efficiency
    // and to avoid round-off error-related problems
    if (nPoints == 3)
    {
        subFaceCentre =
                (1.0 / 3.0) * (subFacePoints[0] + subFacePoints[1] + subFacePoints[2]);

        subFaceArea = 0.5 * ((subFacePoints[1] - subFacePoints[0]) ^
                             (subFacePoints[2] - subFacePoints[0]));
    }
    else if (nPoints > 0)
    {
        vector sumN{Zero};
        scalar sumA{0};
        vector sumAc{Zero};

        point fCentre = subFacePoints[0];
        // initial guess of centre as average of subFacePoints
        for (label pi = 1; pi < nPoints; pi++)
        {
            fCentre += subFacePoints[pi];
        }

        fCentre /= nPoints;

        // loop sub triangles
        for (label pi = 0; pi < nPoints; pi++)
        {
            const point& nextPoint = subFacePoints[(pi + 1) % nPoints];

            vector c = subFacePoints[pi] + nextPoint + fCentre;
            vector n =
                (nextPoint - subFacePoints[pi]) ^ (fCentre - subFacePoints[pi]);
            scalar a = mag(n);

            sumN += n;
            sumA += a;
            sumAc += a * c;
        }

        // This is to deal with zero-area faces. Mark very small faces
        // to be detected in e.g., processorPolyPatch.
        if (sumA < ROOTVSMALL)
        {
            subFaceCentre = fCentre;
            subFaceArea = Zero;
        }
        else
        {
            subFaceCentre = (1.0 / 3.0) * sumAc / sumA;
            subFaceArea = 0.5 * sumN;
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cutFace::cutFace(const fvMesh& mesh)
:
    mesh_(mesh)
{}


// ************************************************************************* //
