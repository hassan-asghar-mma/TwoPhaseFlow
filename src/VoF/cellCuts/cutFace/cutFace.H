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

Class
    Foam::cutFace

Description
    Base class for cutting a face, faceI, of an fvMesh, mesh_, at its
    intersections

    Original code supplied by
        Johan Roenby, DHI (2016)
        Henning Scheufler, DLR (2019)

SourceFiles
    cutFace.C

\*---------------------------------------------------------------------------*/

#ifndef cutFace_H
#define cutFace_H

#include "fvMesh.H"
#include "surfaceFields.H"
#include "volFields.H"
#include <unordered_map>
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class cutFaces Declaration
\*---------------------------------------------------------------------------*/

class cutFace
{
    // Private Data

        //- Reference to mesh
        const fvMesh& mesh_;


protected:

    //- Calculate cut points along edges of face with pointStatus, pointfield
    //- and computes geometric information and the face status where:
    // -  -1: faceStatus is fully below the isosurface
    // -   0: faceStatus is cut (has values larger and smaller than isoValue)
    // -  +1: faceStatus is fully above the isosurface
    void calcSubFace
    (
        const label faceI,
        const scalarList& pointStatus,
        label firstFullySubmergedPoint,
        DynamicList<point>& subFacePoints,
        DynamicList<point>& surfacePoints,
        label& faceStatus,
        vector& subFaceCentre,
        vector& subFaceArea
     );

    //- Calculate cut points along edges of face with pointStatus, pointfield,
    //- weights and computes geometric information and the face status where:
    // -  -1: faceStatus is fully below the isosurface
    // -   0: faceStatus is cut (has values larger and smaller than isoValue)
    // -  +1: faceStatus is fully above the isosurface
    void calcSubFace
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
    );

    //- Calculates centre and normal of the face
    void calcSubFaceCentreAndArea
    (
        DynamicList<point>& subFacePoints,
        vector& subFaceCentre,
        vector& subFaceArea
    );

    //- Calculate cut points along edges of face with pointStatus, pointfield
    //- and computes geometric information and the face status where:
    // -  -1: faceStatus is fully below the isosurface
    // -   0: faceStatus is cut (has values larger and smaller than isoValue)
    // -  +1: faceStatus is fully above the isosurface
    void calcSubFace
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
    );


public:

    // Static Data
    static int debug;

    static int timeIndex_;
    
    // Static map: Key is face ID, value is a list of intersection points
   //static std::unordered_map<std::string, std::unordered_map<Foam::label, std::unordered_map<Foam::label, Foam::point>>> faceIntersectionsMap;
   //static std::unordered_map<std::string, std::unordered_map<Foam::label, std::map<std::pair<Foam::point, Foam::point>, Foam::point>>> faceIntersectionsMap;
   
   static std::unordered_map<Foam::label, std::vector<std::pair<Foam::point, Foam::point>>> patchEdges;
   bool isBottomAdjacent(const Foam::point& edgeStart, const Foam::point& edgeEnd);
    // Function to calculate and store intersections
    static void storeAndCalculateIntersections(
        const fvMesh& mesh_,
        label faceI,
        const scalarList& pointStatus,
        const pointField& points,
        const face& f,
        DynamicList<point>& subFacePoints,
        DynamicList<point>& surfacePoints,
        label firstFullySubmergedPoint
    );

    static bool faceInMap
    (
        label faceI,
        std::unordered_map<Foam::label, std::vector<Foam::point>> faceIntersectionsMap
    );

    void printPatchEdges(const std::unordered_map<label, std::vector<std::pair<point, point>>>& patchEdges) const;
    bool isEdgeInPatchEdges(const std::vector<point>& edge,const std::unordered_map<label, std::vector<std::pair<point, point>>>& patchEdges) const;

    void extractPatchEdges
    (
        const std::string& patchName,
        std::unordered_map<label, std::vector<std::pair<point, point>>>& patchEdges
    );
            
    static std::vector<double> pOld;
    static std::vector<std::vector<point>> faceIntersections; // edge vertices , intersection point makes a vector and the outer vector is for all edges of patch bottom 

    // Constructors

        //- Construct from fvMesh
        explicit cutFace(const fvMesh& mesh);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
