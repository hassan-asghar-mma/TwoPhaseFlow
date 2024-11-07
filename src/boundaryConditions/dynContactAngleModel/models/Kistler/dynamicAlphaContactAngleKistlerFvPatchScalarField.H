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
    Foam::dynamicAlphaContactAngleKistlerFvPatchScalarField

Description:
    This boundary condition employs the Kistler dynamic contact angle model (taken from Dirk Drunding and 
    Anja Lippert dissertation)

    // theta_d = f_Hoff (Ca + f^(-1)_Hoff (\theta0_))
    // f_Hoff (x) := arccos (1- 2*tanh(5.16(x / (1 + 1.31*x^(0.99)))^0.706))

Developed by:
    Muhammad Hassan Asghar
    Mathematical Modeling and Analysis
    TU  Darmstadt

SourceFiles
    dynamicAlphaContactAngleKistlerFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef dynamicAlphaContactAngleKistlerFvPatchScalarField_H
#define dynamicAlphaContactAngleKistlerFvPatchScalarField_H

#include "alphaContactAngleTwoPhaseFvPatchScalarField.H"
#include "volFields.H"
#include "fvMesh.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
               Class Kistler Declaration
\*---------------------------------------------------------------------------*/

class dynamicAlphaContactAngleKistlerFvPatchScalarField
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
    TypeName("dynamicAlphaContactAngleKistler");


    // Constructors

        //- Construct from patch and internal field
        dynamicAlphaContactAngleKistlerFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        dynamicAlphaContactAngleKistlerFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  Kistler
        //  onto a new patch
        dynamicAlphaContactAngleKistlerFvPatchScalarField
        (
            const dynamicAlphaContactAngleKistlerFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        dynamicAlphaContactAngleKistlerFvPatchScalarField
        (
            const dynamicAlphaContactAngleKistlerFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new dynamicAlphaContactAngleKistlerFvPatchScalarField(*this)
            );
        }
        
         //- Evaluate and return dynamic contact-angle
        virtual tmp<scalarField> theta
        (
            const fvPatchVectorField& Up,
            const fvsPatchVectorField& nHat
        ) const;

        //- Construct as copy setting internal field reference
        dynamicAlphaContactAngleKistlerFvPatchScalarField
        (
            const dynamicAlphaContactAngleKistlerFvPatchScalarField&,
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
                new dynamicAlphaContactAngleKistlerFvPatchScalarField(*this, iF)
            );
        }


     //Employs the newton method to solve for the Inverse Hoffman function
       double newton_method
       (
            const std::function<double(double)>& func, 
            const std::function<double(double)>& deriv, 
            scalar initial_guess
       )const;

        // Inverse Hoffman function: Argument in the DCA model 
       double hoffman_function_inv 
       (
            scalar theta_e
       )const;

        // Evaluate the dynamic contact angle based on equilibrium contact angle theta_e
        // and Capillary number at the wall.
       double kistler_model
       (
            scalar theta_e, 
            scalar Ca
       ) const;


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
