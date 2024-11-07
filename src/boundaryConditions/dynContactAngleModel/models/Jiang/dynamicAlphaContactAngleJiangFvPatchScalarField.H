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

Class
    Foam::dynamicAlphaContactAngleJiangFvPatchScalarField

Description
    A dynamic alphaContactAngle scalar boundary condition
    employs the empirical Jiang model for dynamic contact angle 
    ((taken from Dirk Drunding and Anja Lippert dissertation)) i.e.,
    
    // theta = arccos (cos\theta_e - tanh(4.96*Ca^(0.702))(cos\theta_e + 1)) 

Developed by:
    Muhammad Hassan Asghar
    Mathematical Modeling and Analysis
    TU  Darmstadt


SourceFiles
    dynamicAlphaContactAngleJiangFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef dynamicAlphaContactAngleJiangFvPatchScalarField_H
#define dynamicAlphaContactAngleJiangFvPatchScalarField_H

#include "alphaContactAngleTwoPhaseFvPatchScalarField.H"
#include "volFields.H"
#include "fvMesh.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
               Class coxVoinov Declaration
\*---------------------------------------------------------------------------*/

class dynamicAlphaContactAngleJiangFvPatchScalarField
:
    public alphaContactAngleTwoPhaseFvPatchScalarField
{
    // Private data

        //- Equilibrium contact angle
        scalar theta0_;

	//- Visualization 
        mutable volScalarField contactLineAngle_;

    // Private member functions
       
        //- Determine whether a face contains a contact line
        bool hasContactLine
        (
            label faceI
        ) const;


public:

    //- Runtime type information
    TypeName("dynamicAlphaContactAngleJiang");


    // Constructors

        //- Construct from patch and internal field
        dynamicAlphaContactAngleJiangFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        dynamicAlphaContactAngleJiangFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  coxVoinov
        //  onto a new patch
        dynamicAlphaContactAngleJiangFvPatchScalarField
        (
            const dynamicAlphaContactAngleJiangFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        dynamicAlphaContactAngleJiangFvPatchScalarField
        (
            const dynamicAlphaContactAngleJiangFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new dynamicAlphaContactAngleJiangFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        dynamicAlphaContactAngleJiangFvPatchScalarField
        (
            const dynamicAlphaContactAngleJiangFvPatchScalarField&,
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
                new dynamicAlphaContactAngleJiangFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        //- Evaluate and return dynamic contact-angle
        virtual tmp<scalarField> theta
        (
            const fvPatchVectorField& Up,
            const fvsPatchVectorField& nHat
        ) const;


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
