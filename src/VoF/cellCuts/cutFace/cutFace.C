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
//std::unordered_map<Foam::label, std::vector<Foam::point>> Foam::cutFace::faceIntersectionsMap;
std::unordered_map<std::string, std::unordered_map<Foam::label, std::unordered_map<Foam::label, Foam::point>>> Foam::cutFace::faceIntersectionsMap;
std::unordered_map<Foam::label, std::vector<std::pair<Foam::point, Foam::point>>> Foam::cutFace::patchEdges;

int Foam::cutFace::timeIndex_ = -100;
// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * *

bool Foam::cutFace::faceInMap(label faceI,std::unordered_map<Foam::label, std::vector<Foam::point>> faceIntersectionsMap)
{   
    return faceIntersectionsMap.find(faceI) != faceIntersectionsMap.end();
}

//Extracting all the edges of the patch
//Needed to remove the repition during saving of all intersection points in bottom patch
// Idea is to check if the edge from other patches coincidde with edge of this container, do not calculate the intersection point
// rather just retrieve the priously calculated intersection point : avind repition and make sure that other patches also have the same intersection points for bottom faces
// is used in storeAndCalculateIntersections function
void Foam::cutFace::extractPatchEdges
(
    const std::string& patchName,
    std::unordered_map<label, std::vector<std::pair<point, point>>>& patchEdges
)
{
    // Check if the map is already filled
    if (!patchEdges.empty())
    {
        return; // Skip if the edges are already extracted
    }

    // Locate the patch by name
    label patchID = -1;
    for (label i = 0; i < mesh_.boundaryMesh().size(); ++i)
    {
        if (mesh_.boundaryMesh()[i].name() == patchName)
        {
            patchID = i;
            break;
        }
    }

    // Access the specified patch
    const auto& patch = mesh_.boundaryMesh()[patchID];
    const pointField& points = mesh_.points();

    // Iterate through faces of the patch
    for (label localFaceID = 0; localFaceID < patch.size(); ++localFaceID)
    {
        // Get global face ID
        label globalFaceID = patch.start() + localFaceID;

        const face& patchFace = mesh_.faces()[globalFaceID];

        // Debug: Output points of the face
        Info << "Global Face ID: " << globalFaceID << " Points: ";
        for (label i = 0; i < patchFace.size(); ++i)
        {
            const point& p = points[patchFace[i]];
            Info << "(" << p.x() << ", " << p.y() << ", " << p.z() << ") ";
        }
        Info << endl;

        // Construct edges for the face
        for (label i = 0; i < patchFace.size(); ++i)
        {
            label startIdx = patchFace[i];                       // Global start index
            label endIdx = patchFace[(i + 1) % patchFace.size()]; // Global end index

            const point& startPoint = points[startIdx];
            const point& endPoint = points[endIdx];

            // Debug: Output each edge
            Info << "Edge: (" << startPoint.x() << ", " << startPoint.y() << ", " << startPoint.z() << ") -> "
                 << "(" << endPoint.x() << ", " << endPoint.y() << ", " << endPoint.z() << ")" << endl;

            // Add edge to the map using the global face ID
            patchEdges[globalFaceID].emplace_back(startPoint, endPoint);
        }
    }
}


void Foam::cutFace::storeAndCalculateIntersections(
    const fvMesh& mesh_,
    label faceI,
    const scalarList& pointStatus,
    const pointField& points,
    const face& f,
    DynamicList<point>& subFacePoints,
    DynamicList<point>& surfacePoints,
    label firstFullySubmergedPoint
)
{

    label patchi = 0;
    label bottomPatchId = 0;

    for (label patchID = 0; patchID < mesh_.boundaryMesh().size(); ++patchID)
    {
        const auto patch = mesh_.boundaryMesh()[patchID];
        const label patchStart = patch.start();
        const label patchEnd = patchStart + patch.size();
        if (mesh_.boundaryMesh()[patchID].name()=="bottom")
        {
            bottomPatchId = patchID;
        }

        // Check if the face ID falls within this patch's range
        if (faceI >= patchStart && faceI < patchEnd)
        {
            patchi = patchID;
        }
    }

    const polyBoundaryMesh& boundaryMesh = mesh_.boundaryMesh();
    const label start = boundaryMesh[patchi].start();
    word patchName = mesh_.boundaryMesh()[patchi].name();

    

    Info << "#########################################\n";
    // Info << "For face " << faceI << " with local face id " << faceI - start 
    //         <<  " in the patch " << patchName << " I will be running for " << f.size()- firstFullySubmergedPoint << " times \n";

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
            Info << "\n" << i << " ith iteration: \n cut point status in first else if is the interface is on vertex " << points[f[idx]] << endl;
            subFacePoints.append(points[f[idx]]);
            surfacePoints.append(points[f[idx]]);
        }

        // Check if this edge coincides with a bottom patch edge
        bool isBottomAdjacent = false;

        // Skip saving intersection points for non-bottom patches unless edges coincide with the bottom patch
        if (patchName != "bottom")
        {
            for (const auto& bottomFaceEntry : patchEdges)
            {
                for (const auto& edge : bottomFaceEntry.second)
                {
                    // Compare edges (direct and reversed)
                    if ((edge.first == points[f[idx]] && edge.second == points[f[nextIdx]]) ||
                        (edge.first == points[f[nextIdx]] && edge.second == points[f[idx]]))
                    {
                        isBottomAdjacent = true;
                        break;
                    }
                }
                if (isBottomAdjacent)
                {
                    break;
                }
            }
        }


        point intersection;

        //If the intersection point is already present in the map
        //Check for the presence of patch, patch face, and then corresponding intersection point
        if (faceIntersectionsMap.find(patchName) != faceIntersectionsMap.end() &&
            faceIntersectionsMap[patchName].find(faceI) != faceIntersectionsMap[patchName].end() &&
            faceIntersectionsMap[patchName][faceI].find(idx) != faceIntersectionsMap[patchName][faceI].end())
        {
            // Use the stored intersection point
            intersection = faceIntersectionsMap[patchName][faceI][idx];
            subFacePoints.append(intersection);
            surfacePoints.append(intersection);
            Info << "Retrieved intersection for Patch: " << patchName
                 << ", Face: " << faceI << ", Edge: " << idx
                 << ": " << intersection << endl;
            continue; // Skip to the next edge
        }   

        // The intersection point is not present in the map and should calculated and stored 
        // Calculate intersection points for edges
        if (
            (pointStatus[idx] < 0 && pointStatus[nextIdx] > 0) ||
            (pointStatus[idx] > 0 && pointStatus[nextIdx] < 0)
            )
        {
            // Compute intersection point on the edge
            label nextP = f.nextLabel(idx);
            vector dir = points[nextP] - points[f[idx]];
            scalar weight = (0.0 - pointStatus[idx]) /
                            (pointStatus[nextIdx] - pointStatus[idx]);

            intersection = points[f[idx]] + weight * dir;
            Info << "Computed intersection for Patch: " << patchName
                 << ", Face: " << faceI << ", Edge: " << idx
                 << ": " << intersection << endl;

            // Store the intersection point in the map if relevant
            if (patchName == "bottom" || isBottomAdjacent)
            {
                faceIntersectionsMap[patchName][faceI][idx] = intersection;
            }


            // Append intersection points to subFacePoints and surfacePoints
            subFacePoints.append(intersection);
            surfacePoints.append(intersection);

            // Save the intersection point to the map
            //Info << "Surface points: " <<surfacePoints  << endl;


        }
    }

    // if(faceIntersectionsMap.empty())
    // {
    //     Info << "The container is empty" << endl;
    //     faceIntersectionsMap[faceI] = intersectionPoints;
    // }
    // Store the intersections in the static map
    if (!faceIntersectionsMap[patchName][faceI].empty())
    {   
        Info << "#########################################\n";
        Info << "Face ID: " << faceI << ", Patch: " << patchName << endl;
        Info << "SurfacePoints are:" << "\n__________________________" << "\n"  << surfacePoints << endl;
        Info << "Intersection Points: ";
        for (const auto& edgePoint : faceIntersectionsMap[patchName][faceI])
        {
            Info << "\nEdge Index " << edgePoint.first << ": "
                 << "(" << edgePoint.second.x() << ", "
                 << edgePoint.second.y() << ", "
                 << edgePoint.second.z() << ")";
                 
        }
        Info << "\n#########################################\n";
    }


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

    //Compute the bottom patch edges
    extractPatchEdges("bottom", patchEdges);

    // Call the function to calculate and store intersections
    storeAndCalculateIntersections(mesh_, faceI, pointStatus, points, f, subFacePoints, surfacePoints, firstFullySubmergedPoint);

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
