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
//std::unordered_map<std::vector<Foam::point>> Foam::cutFace::faceIntersectionsMap;
//std::unordered_map<std::string, std::unordered_map<Foam::label, std::unordered_map<Foam::label, Foam::point>>> Foam::cutFace::faceIntersectionsMap;

std::unordered_map<Foam::label, std::vector<std::pair<Foam::point, Foam::point>>> Foam::cutFace::patchEdges; //store all the edges before, 
int Foam::cutFace::timeIndex_ = -100;
std::vector<std::vector<Foam::point>> Foam::cutFace::faceIntersections;

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * *

bool Foam::cutFace::faceInMap(label faceI,std::unordered_map<Foam::label, std::vector<Foam::point>> faceIntersectionsMap)
{   
    return faceIntersectionsMap.find(faceI) != faceIntersectionsMap.end();
}

void Foam::cutFace::printPatchEdges(const std::unordered_map<label, std::vector<std::pair<point, point>>>& patchEdges) const
{
    Info << "Printing edges of all patches:" << endl;
    for (const auto& patchEntry : patchEdges)
    {
        label faceID = patchEntry.first;
        const auto& edges = patchEntry.second;

        Info << "Face local id: " << faceID << endl;
        for (const auto& edge : edges)
        {
            const point& start = edge.first;
            const point& end = edge.second;

            Info << "  Edge: (" << start.x() << ", " << start.y() << ", " << start.z() << ") -> "
                 << "(" << end.x() << ", " << end.y() << ", " << end.z() << ")" << endl;
        }
    }
    Info << "Finished printing edges of all patches." << endl;
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
        // Info << "Global Face ID: " << globalFaceID << " Points: ";
        // for (label i = 0; i < patchFace.size(); ++i)
        // {
        //     const point& p = points[patchFace[i]];
        //     Info << "(" << p.x() << ", " << p.y() << ", " << p.z() << ") ";
        // }
        // Info << endl;

        // Construct edges for the face
        for (label i = 0; i < patchFace.size(); ++i)
        {
            label startIdx = patchFace[i];                       // Global start index
            label endIdx = patchFace[(i + 1) % patchFace.size()]; // Global end index

            const point& startPoint = points[startIdx];
            const point& endPoint = points[endIdx];

            // Debug: Output each edge
            // Info << "Edge: (" << startPoint.x() << ", " << startPoint.y() << ", " << startPoint.z() << ") -> "
            //      << "(" << endPoint.x() << ", " << endPoint.y() << ", " << endPoint.z() << ")" << endl;

            // Add edge to the map using the global face ID
            patchEdges[localFaceID].emplace_back(startPoint, endPoint);
        }
    }
}

bool Foam::cutFace::isEdgeInPatchEdges(
    const std::vector<point>& edge,
    const std::unordered_map<label, std::vector<std::pair<point, point>>>& patchEdges
) const
{
    // Ensure the edge is ordered consistently (start point is lexicographically smaller than end point)
    point start = edge[0];
    point end = edge[1];
 
    for (const auto& patchEntry : patchEdges)
    {
        const auto& edges = patchEntry.second;
        for (const auto& storedEdge : edges)
        {
            // Compare both directions of the edge
            if ((storedEdge.first == start && storedEdge.second == end) ||
                (storedEdge.first == end && storedEdge.second == start))
            {
                return true; // Edge found
            }
        }
    }
    return false; // Edge not found
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
    //Store the bottom edges at the first instance
    if(mesh_.time().timeIndex() == 0 && patchEdges.empty())
    {
        extractPatchEdges("bottom", patchEdges);
    }
    // for (int count=0; count<1;count++)
    //     printPatchEdges(patchEdges);
    
    const pointField& points = mesh_.points();
    const face& f = mesh_.faces()[faceI];
    const auto & faceOwner = mesh_.faceOwner();
    Info << "\n###################################\n";
    Info << "Local Face ID " << faceI <<" of cell "<< faceOwner[faceI]<< endl;
    if (firstFullySubmergedPoint == -1) // is in gasPhase
    {
        faceStatus = 1;
        Info << "Face " << faceI << " is in gas phase\n";
        subFaceCentre = Zero;
        subFaceArea = Zero;
        return;
    }

    // loop face and append the cuts
    // loop starts at firstFullySubmergedPoint
    static int counter = 0;
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

        //Current Edge
        std::vector<point> edge = {points[f.nextLabel(idx)], points[f[idx]]};
        //Info <<" Edge " << edge << " is present? "<< isEdgeInPatchEdges(edge, patchEdges) << endl;

        if (pointStatus[idx] > 0) // append fluid point
        {
            Info << "Appending vertex " << points[f[idx]] <<" . This is in liquid\n";
            subFacePoints.append(points[f[idx]]);
        }
        else if (pointStatus[idx] == 0) // append cut point
        {
            Info << "Appending vertex " << points[f[idx]] <<" . This is on vertex\n";
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
            Info <<"i: " << i <<  " Face " << faceI << " has intersection on edge: \n\t" << points[nextP] << " -> " << points[f[idx]] << endl;
            scalar weight =
                (0.0 - pointStatus[idx]) /
                (pointStatus[nextIdx] - pointStatus[idx]); // cutValue is zero

            point p = points[f[idx]] + weight * dir;

            if (mesh_.time().timeIndex() == 0 && isEdgeInPatchEdges(edge, patchEdges))
            {
                bool edgeAlreadyExists = false;

                // Iterate through the existing faceIntersections to check for the edge
                for (const auto& intersection : faceIntersections)
                {
                    // Ensure the edge matches in either order
                    if ((intersection[0] == edge[0] && intersection[1] == edge[1]) ||
                        (intersection[0] == edge[1] && intersection[1] == edge[0]))
                    {
                        edgeAlreadyExists = true;
                        break;
                    }
                }

                // Add the new intersection if the edge is not already present
                if (!edgeAlreadyExists)
                {
                    faceIntersections.push_back({edge[0], edge[1], p});
                    counter++;
                }
            }
            if (mesh_.time().timeIndex() >=0 && isEdgeInPatchEdges(edge, patchEdges))
            {
                bool edgeAlreadyExists = false;

                // Iterate through the existing faceIntersections to check for the edge
                for (const auto& intersection : faceIntersections)
                {
                    // Ensure the edge matches in either order
                    if ((intersection[0] == edge[0] && intersection[1] == edge[1]) ||
                        (intersection[0] == edge[1] && intersection[1] == edge[0]))
                    {
                        edgeAlreadyExists = true;
                        p = intersection[2];
                        break;
                    }
                }
            }

        
            Info << "\t\tIntersection point : " << p << endl << endl;
            subFacePoints.append(p);
            surfacePoints.append(p);
        }
    }
    Info << "Subface is: " << subFacePoints<< endl;
    Info << "Surface is: " << surfacePoints << endl;

    if (subFacePoints.size() >= 3)
    {
        faceStatus = 0;
        calcSubFaceCentreAndArea(subFacePoints, subFaceCentre, subFaceArea);
    }
    else
    {
        faceStatus = -1;
    }
    // for (size_t i = 0; i < faceIntersections.size(); ++i) // Iterate over the outer vector
    // {
    //     const auto& innerVector = faceIntersections[i];
    //     Info << "Face " << i << " edge-edge-intersections:" << endl;

    //     for (size_t j = 0; j < innerVector.size(); ++j) // Iterate over the inner vector
    //     {
    //         const auto& point = innerVector[j];
    //         Info << "    Point " << j << ": " << point << endl;
    //     }
    // }

}
// void Foam::cutFace::storeAndCalculateIntersections(
//     const fvMesh& mesh_,
//     label faceI,
//     const scalarList& pointStatus,
//     const pointField& points,
//     const face& f,
//     DynamicList<point>& subFacePoints,
//     DynamicList<point>& surfacePoints,
//     label firstFullySubmergedPoint
// )
// {


//     label patchi = 0;
//     label bottomPatchId = 0;

//     for (label patchID = 0; patchID < mesh_.boundaryMesh().size(); ++patchID)
//     {
//         const auto patch = mesh_.boundaryMesh()[patchID];
//         const label patchStart = patch.start();
//         const label patchEnd = patchStart + patch.size();
//         if (mesh_.boundaryMesh()[patchID].name() == "bottom")
//         {
//             bottomPatchId = patchID;
//         }

//         // Check if the face ID falls within this patch's range
//         if (faceI >= patchStart && faceI < patchEnd)
//         {
//             patchi = patchID;
//         }
//     }

//     const polyBoundaryMesh& boundaryMesh = mesh_.boundaryMesh();
//     const label start = boundaryMesh[patchi].start();
//     word patchName = mesh_.boundaryMesh()[patchi].name();

//     Info << "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&\n";
//     Info << "For face " << faceI << " with local face id " << faceI - start
//          << " in the patch " << patchName << " I am running\n\n";

//     // Check if the time index is zero for initialization phase
//     bool allowUpdate = (mesh_.time().timeIndex() == 0);

//     // Loop starting at firstFullySubmergedPoint
//     for (label i = firstFullySubmergedPoint; i < firstFullySubmergedPoint + f.size(); ++i)
//     {

//         label idx = i % f.size();
//         label nextIdx = (i + 1) % f.size();

//         // Append fluid points
//         if (pointStatus[idx] > 0)
//         {
//             subFacePoints.append(points[f[idx]]);
//         }
//         // Append cut points
//         else if (pointStatus[idx] == 0)
//         {
//             subFacePoints.append(points[f[idx]]);
//             surfacePoints.append(points[f[idx]]);
//         }

//         // Check if this edge coincides with a bottom patch edge
//         bool isBottomAdjacent = false;
//         if (patchName != "bottom" || mesh_.isInternalFace(faceI))
//         {
//             for (const auto& bottomFaceEntry : patchEdges)
//             {
//                 for (const auto& edge : bottomFaceEntry.second)
//                 {
//                     if ((edge.first == points[f[idx]] && edge.second == points[f[nextIdx]]) ||
//                         (edge.first == points[f[nextIdx]] && edge.second == points[f[idx]]))
//                     {
//                         isBottomAdjacent = true;
//                         Info << "The above Edge is adjacent to the bottom\n";
//                         break;
//                     }
//                 }
//                 if (isBottomAdjacent)
//                 {
//                     break;
//                 }
//             }
//         }

//         point intersection;

//         // If the intersection point is already present in the map
//         if (faceIntersectionsMap.find(patchName) != faceIntersectionsMap.end() &&
//             faceIntersectionsMap[patchName].find(faceI) != faceIntersectionsMap[patchName].end() &&
//             faceIntersectionsMap[patchName][faceI].find(idx) != faceIntersectionsMap[patchName][faceI].end())
//         {
//             intersection = faceIntersectionsMap[patchName][faceI][idx];
//             subFacePoints.append(intersection);
//             surfacePoints.append(intersection);
//             Info << "Retrieved point " << intersection << " from " << patchName << " - " << faceI << " - " << idx << endl;
//             continue; // Skip to the next edge
//         }

//         // Calculate intersection points for edges
//         if ((pointStatus[idx] < 0 && pointStatus[nextIdx] > 0) ||
//             (pointStatus[idx] > 0 && pointStatus[nextIdx] < 0))
//         {
//             // Compute intersection point on the edge
//             label nextP = f.nextLabel(idx);
//             vector dir = points[nextP] - points[f[idx]];
//             scalar weight = (0.0 - pointStatus[idx]) / (pointStatus[nextIdx] - pointStatus[idx]);
//             intersection = points[f[idx]] + weight * dir;

//             // Store the intersection point in the map if relevant and updates are allowed
//             if (allowUpdate && (patchName == "bottom" || isBottomAdjacent))
//             {
//                 Info << "Saving point " << intersection << " in " << patchName << " - " << faceI << " - " << idx << endl;
//                 faceIntersectionsMap[patchName][faceI][idx] = intersection;
//             }

//             subFacePoints.append(intersection);
//             surfacePoints.append(intersection);
//         }
//     }

//     Info << "END&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&\n";

//     if (!faceIntersectionsMap[patchName][faceI].empty())
//     {
//         Info << "\n\n\n#########################################\n";
//         Info << "Face ID: " << faceI << ", Patch: " << patchName << endl;
//         Info << "SurfacePoints are:\n__________________________\n" << surfacePoints << endl;
//         Info << "Intersection Points: ";
//         for (const auto& edgePoint : faceIntersectionsMap[patchName][faceI])
//         {
//             Info << "\nEdge Index " << edgePoint.first << ": "
//                  << "(" << edgePoint.second.x() << ", "
//                  << edgePoint.second.y() << ", "
//                  << edgePoint.second.z() << ")";
//         }
//         Info << "\n#########################################\n";
//     }
// }


// void Foam::cutFace::calcSubFace
// (
//     const label faceI,
//     const scalarList& pointStatus,
//     label firstFullySubmergedPoint,
//     DynamicList<point>& subFacePoints,
//     DynamicList<point>& surfacePoints,
//     label& faceStatus,
//     vector& subFaceCentre,
//     vector& subFaceArea
// )
// {
//     const pointField& points = mesh_.points();
//     const face& f = mesh_.faces()[faceI];

//     if (firstFullySubmergedPoint == -1) // In gas phase
//     {
//         faceStatus = 1;
//         subFaceCentre = Zero;
//         subFaceArea = Zero;
//         return;
//     }

//     //Compute the bottom patch edges
//     extractPatchEdges("bottom", patchEdges);
//     Info << "Extraction of bottom patch faces is successfull\n\n";

//     Info << timeIndex_ << " ---  " <<  mesh_.time().timeIndex()<< "\n";
//     // Call the function to calculate and store intersections
//     storeAndCalculateIntersections(mesh_, faceI, pointStatus, points, f, subFacePoints, surfacePoints, firstFullySubmergedPoint);

//     // Additional logic for sub-face area and center
//     if (subFacePoints.size() >= 3)
//     {
//         faceStatus = 0;
//         calcSubFaceCentreAndArea(subFacePoints, subFaceCentre, subFaceArea);
//     }
//     else
//     {
//         faceStatus = -1;
//     }
// }
// void Foam::cutFace::calcSubFace(
//     const label faceI,
//     const scalarList& pointStatus,
//     label firstFullySubmergedPoint,
//     DynamicList<point>& subFacePoints,
//     DynamicList<point>& surfacePoints,
//     label& faceStatus,
//     vector& subFaceCentre,
//     vector& subFaceArea
// )
// {
//     const pointField& points = mesh_.points();
//     const face& f = mesh_.faces()[faceI];

//     if (firstFullySubmergedPoint == -1) // is in gasPhase
//     {
//         faceStatus = 1;
//         subFaceCentre = Zero;
//         subFaceArea = Zero;
//         return;
//     }

//     bool allowUpdate = (mesh_.time().timeIndex() == 0);

//     // Determine the patch name and face ID
//     label patchi = 0;
//     for (label patchID = 0; patchID < mesh_.boundaryMesh().size(); ++patchID)
//     {
//         const auto patch = mesh_.boundaryMesh()[patchID];
//         const label patchStart = patch.start();
//         const label patchEnd = patchStart + patch.size();
//         if (faceI >= patchStart && faceI < patchEnd)
//         {
//             patchi = patchID;
//             break;
//         }
//     }

//     const word patchName = mesh_.boundaryMesh()[patchi].name();

//     // Loop through face edges
//     for (label i = firstFullySubmergedPoint; i < firstFullySubmergedPoint + f.size(); ++i)
//     {
//         label idx = i % f.size();
//         label nextIdx = (i + 1) % f.size();

//         // Normalize the edge for consistent representation
//         point edgeStart = points[f[idx]];
//         point edgeEnd = points[f[nextIdx]];
//         if (edgeStart > edgeEnd)
//         {
//             std::swap(edgeStart, edgeEnd);
//         }

//         Info << "Processing face " << faceI << ", edge: "
//              << edgeStart << " -> " << edgeEnd << endl;

//         std::pair<point, point> edgeKey = {edgeStart, edgeEnd};
//         point intersection;

//         // Check if the edge exists in the map for the current patch or the bottom patch
//         for (const auto& patch : std::vector<std::string>{patchName, "bottom"})
//         {
//             auto patchEntry = faceIntersectionsMap.find(patch);
//             if (patchEntry != faceIntersectionsMap.end())
//             {
//                 auto faceEntry = patchEntry->second.find(faceI);
//                 if (faceEntry != patchEntry->second.end())
//                 {
//                     for (const auto& edgeData : faceEntry->second)
//                     {
//                         // Check if the edge matches the edge under observation
//                         if ((edgeData.first.first == edgeStart && edgeData.first.second == edgeEnd) ||
//                             (edgeData.first.first == edgeEnd && edgeData.first.second == edgeStart))
//                         {
//                             intersection = edgeData.second;
//                             Info << "Using stored intersection for face " << faceI
//                                  << " edge " << idx << ": " << intersection << endl;
//                             subFacePoints.append(intersection);
//                             surfacePoints.append(intersection);
//                             goto nextEdge;
//                         }
//                     }
//                 }
//             }
//         }

//         // Calculate intersection if not found in the map
//         if (
//             (pointStatus[idx] < 0 && pointStatus[nextIdx] > 0) ||
//             (pointStatus[idx] > 0 && pointStatus[nextIdx] < 0))
//         {
//             vector dir = edgeEnd - edgeStart;
//             scalar weight = (0.0 - pointStatus[idx]) /
//                             (pointStatus[nextIdx] - pointStatus[idx]);
//             intersection = edgeStart + weight * dir;
//             Info << "Calculated intersection for edge: " 
//                  << edgeStart << " -> " << edgeEnd 
//                  << " is: " << intersection << endl;

//             if (allowUpdate && (patchName == "bottom" || isBottomAdjacent(edgeStart, edgeEnd)))
//             {
//                 faceIntersectionsMap[patchName][faceI][edgeKey] = intersection;
//                 Info << "Storing intersection for face " << faceI
//                      << " edge " << idx << ": " << intersection << endl;
//             }

//             subFacePoints.append(intersection);
//             surfacePoints.append(intersection);
//         }

//         Info << "Stored intersections for face " << faceI << ": " << endl;
//         for (const auto& edgeData : faceIntersectionsMap[patchName][faceI])
//         {
//             Info << "Edge: " << edgeData.first.first << " -> " << edgeData.first.second
//                  << ", Intersection: " << edgeData.second << endl;
//         }

//     nextEdge:
//         continue;
//     }

//     if (subFacePoints.size() >= 3)
//     {
//         faceStatus = 0;
//         calcSubFaceCentreAndArea(subFacePoints, subFaceCentre, subFaceArea);
//     }
//     else
//     {
//         faceStatus = -1;
//     }

//     Info << "\n__________________________";
//     Info << "Face " << faceI << " SurfacePoints are:" << "\n__________________________" << "\n" << surfacePoints << endl;
// }

// bool Foam::cutFace::isBottomAdjacent(const point& edgeStart, const point& edgeEnd)
// {
//     for (const auto& bottomFaceEntry : patchEdges)
//     {
//         for (const auto& edge : bottomFaceEntry.second)
//         {
//             if ((edge.first == edgeStart && edge.second == edgeEnd) ||
//                 (edge.first == edgeEnd && edge.second == edgeStart))
//             {
//                 return true;
//             }
//         }
//     }
//     return false;
// }


// void Foam::cutFace::calcSubFace
// (
//     const label faceI,
//     const scalarList& pointStatus,
//     label firstFullySubmergedPoint,
//     DynamicList<point>& subFacePoints,
//     DynamicList<point>& surfacePoints,
//     label& faceStatus,
//     vector& subFaceCentre,
//     vector& subFaceArea
// )
// {
//     const pointField& points = mesh_.points();
//     const face& f = mesh_.faces()[faceI];
//     const auto & faceOwner = mesh_.faceOwner();
//     Info << "\n###################################\n";
//     Info << "Local Face ID " << faceI <<" of cell "<< faceOwner[faceI]<< endl;
//     if (firstFullySubmergedPoint == -1) // is in gasPhase
//     {
//         faceStatus = 1;
//         Info << "Face " << faceI << " is in gas phase\n";
//         subFaceCentre = Zero;
//         subFaceArea = Zero;
//         return;
//     }

//     // loop face and append the cuts
//     // loop starts at firstFullySubmergedPoint
//     for
//     (
//         label i = firstFullySubmergedPoint;
//         i < firstFullySubmergedPoint + f.size();
//         ++i
//     )
//     {
//         // max two points are appended during one cycle
//         label idx = i % f.size();
//         label nextIdx = (i + 1) % f.size();

//         if (pointStatus[idx] > 0) // append fluid point
//         {
//             Info << "Appending vertex " << points[f[idx]] <<" . This is in liquid\n";
//             subFacePoints.append(points[f[idx]]);
//         }
//         else if (pointStatus[idx] == 0) // append cut point
//         {
//             Info << "Appending vertex " << points[f[idx]] <<" . This is on vertex\n";
//             subFacePoints.append(points[f[idx]]);
//             surfacePoints.append(points[f[idx]]);
//         }

//         if
//         (
//             (pointStatus[idx] < 0 && pointStatus[nextIdx] > 0) ||
//             (pointStatus[idx] > 0 && pointStatus[nextIdx] < 0)
//         ) // cut on edge cut Value is zero
//         {
//             label nextP = f.nextLabel(idx);
//             vector dir = points[nextP] - points[f[idx]];
//             Info <<"i: " << i <<  " Face " << faceI << " has intersection on edge: \n\t" << points[nextP] << " -> " << points[f[idx]] << endl;
//             scalar weight =
//                 (0.0 - pointStatus[idx]) /
//                 (pointStatus[nextIdx] - pointStatus[idx]); // cutValue is zero

//             point p = points[f[idx]] + weight * dir;
//             Info << "\t\tIntersection point : " << p << endl << endl;
//             subFacePoints.append(p);
//             surfacePoints.append(p);
//         }
//     }
//     Info << "Subface is: " << subFacePoints<< endl;
//     Info << "Surface is: " << surfacePoints << endl;

//     if (subFacePoints.size() >= 3)
//     {
//         faceStatus = 0;
//         calcSubFaceCentreAndArea(subFacePoints, subFaceCentre, subFaceArea);
//     }
//     else
//     {
//         faceStatus = -1;
//     }
// }



// void Foam::cutFace::calcSubFace
// (
//     const label faceI,
//     const scalarList& pointStatus,
//     label firstFullySubmergedPoint,
//     DynamicList<point>& subFacePoints,
//     DynamicList<point>& surfacePoints,
//     label& faceStatus,
//     vector& subFaceCentre,
//     vector& subFaceArea
// )
// {
//     const pointField& points = mesh_.points();
//     const face& f = mesh_.faces()[faceI];

//     if (firstFullySubmergedPoint == -1) // is in gasPhase
//     {
//         faceStatus = 1;
//         subFaceCentre = Zero;
//         subFaceArea = Zero;
//         return;
//     }

//     // loop face and append the cuts
//     // loop starts at firstFullySubmergedPoint
//     for
//     (
//         label i = firstFullySubmergedPoint;
//         i < firstFullySubmergedPoint + f.size();
//         ++i
//     )
//     {
//         // max two points are appended during one cycle
//         label idx = i % f.size();
//         label nextIdx = (i + 1) % f.size();
//         if(!mesh_.isInternalFace(faceI))
//             Info << "The current for Edge is: " <<  points[f.nextLabel(idx)] <<" -> " << points[f[idx]] << endl;
//         if (pointStatus[idx] > 0) // append fluid point
//         {
//                 Info << "idx>0\n";      
//             subFacePoints.append(points[f[idx]]);
//         }
//         else if (pointStatus[idx] == 0) // append cut point
//         {
//                 Info << "idx status ==0\n";
//             subFacePoints.append(points[f[idx]]);
//             surfacePoints.append(points[f[idx]]);
//         }

//         if
//         (
//             (pointStatus[idx] < 0 && pointStatus[nextIdx] > 0) ||
//             (pointStatus[idx] > 0 && pointStatus[nextIdx] < 0)
//         ) // cut on edge cut Value is zero
//         {
//             label nextP = f.nextLabel(idx);
//             vector dir = points[nextP] - points[f[idx]];
//             scalar weight =
//                 (0.0 - pointStatus[idx]) /
//                 (pointStatus[nextIdx] - pointStatus[idx]); // cutValue is zero

//             point p = points[f[idx]] + weight * dir;

//             Info << "Point p << " <<p << " in " << points[nextP] << " -> " <<  points[f[idx]] << " for face " << faceI << endl;
//             subFacePoints.append(p);
//             surfacePoints.append(p);
//         }
//     }

//     if (subFacePoints.size() >= 3)
//     {
//         faceStatus = 0;
//         calcSubFaceCentreAndArea(subFacePoints, subFaceCentre, subFaceArea);
//     }
//     else
//     {
//         faceStatus = -1;
//     }
//     Info << "\n__________________________" ;
//     Info << "Face " << faceI << " SurfacePoints are:" << "\n__________________________" << "\n"  << surfacePoints << endl;
// }

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
