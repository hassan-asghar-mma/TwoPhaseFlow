/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2021 OpenCFD Ltd.
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
    Foam::navierSlipHysteresis

Group
    grpWallBoundaryConditions grpGenericBoundaryConditions

Description
    This boundary condition provides a partial slip condition.  The amount of
    slip is controlled by a user-supplied field.

Usage
    \table
      Property      | Description                     | Reqd | Default
      type          | Type name: partialSlip          | yes  | -
      refValue      | Reference value at zero slip    | no   | 0
      valueFraction | Fraction of refValue used for boundary [0-1] | yes |
      writeValue    | Output patch value (eg, ParaView) | no | false
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            partialSlip;
        refValue        uniform 0.001;
        valueFraction   uniform 0.1;
        value           uniform 0;
    }
    \endverbatim

See also
    Foam::transformFvPatchField

SourceFiles
    navierSlipHysteresis.C

\*---------------------------------------------------------------------------*/

#ifndef navierSlipHysteresis_H
#define navierSlipHysteresis_H

#include "transformFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class navierSlipHysteresis Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class navierSlipHysteresis
:
    public transformFvPatchField<Type>
{
    //- The parent boundary condition type
    typedef transformFvPatchField<Type> parent_bctype;


    // Private Data

        //- Reference value
        Field<Type> refValue_;

        //- Fraction (0-1) of value used for boundary condition
        scalarField valueFraction_;

        //- Flag to output patch values (e.g. for ParaView)
        bool writeValue_;

        // threshold for pinning
        scalar gamma_;


public:

    //- Runtime type information
    TypeName("partialSlipHysteresis");


    // Constructors

        //- Construct from patch and internal field
        navierSlipHysteresis
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        navierSlipHysteresis
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given navierSlipHysteresis onto a new patch
        navierSlipHysteresis
        (
            const navierSlipHysteresis<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        navierSlipHysteresis
        (
            const navierSlipHysteresis<Type>&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<Type>> clone() const
        {
            return tmp<fvPatchField<Type>>
            (
                new navierSlipHysteresis<Type>(*this)
            );
        }

        //- Construct as copy setting internal field reference
        navierSlipHysteresis
        (
            const navierSlipHysteresis<Type>&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<Type>> clone
        (
            const DimensionedField<Type, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<Type>>
            (
                new navierSlipHysteresis<Type>(*this, iF)
            );
        }


    // Member Functions

        //- False: this patch field is not altered by assignment
        virtual bool assignable() const { return false; }

        // Access

            virtual Field<Type>& refValue()
            {
                return refValue_;
            }

            virtual const Field<Type>& refValue() const
            {
                return refValue_;
            }

            virtual scalarField& valueFraction()
            {
                return valueFraction_;
            }

            virtual const scalarField& valueFraction() const
            {
                return valueFraction_;
            }

            virtual scalar& gamma()
            {
                return gamma_;
            }

            virtual const scalar& gamma() const
            {
                return gamma_;
            }

        // Function to compute f(x) for viscous stress
            Foam::scalarField fStress(
                const Foam::vectorField& tangentialStress, // stress values for the patch
                const Foam::dimensionedScalar gamma, // limit for slip
                const Foam::scalar f_max, // maximum value of the function = 1.0
                const Foam::dimensionedScalar x_s, // the value of Stress where the function f=f_max; x_s = gamma +\delta
                const Foam::scalar n // exponent value to control the saturation position, must be b/w 1 and 2 to get monotonicity.
            );

            Foam::scalarField computeValueFraction() const;



        // Mapping

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchField<Type>&,
                const labelList&
            );


        // Evaluation

            //- Return gradient at boundary
            virtual tmp<Field<Type>> snGrad() const;

            //- Evaluate the patch field
            virtual void evaluate
            (
                const Pstream::commsTypes commsType =
                    Pstream::commsTypes::blocking
            );

            //- Return face-gradient transform diagonal
            virtual tmp<Field<Type>> snGradTransformDiag() const;







        //- Write
        virtual void write(Ostream&) const;


    // Member Operators

        virtual void operator=(const UList<Type>&) {}

        virtual void operator=(const fvPatchField<Type>&) {}
        virtual void operator+=(const fvPatchField<Type>&) {}
        virtual void operator-=(const fvPatchField<Type>&) {}
        virtual void operator*=(const fvPatchField<scalar>&) {}
        virtual void operator/=(const fvPatchField<scalar>&) {}

        virtual void operator+=(const Field<Type>&) {}
        virtual void operator-=(const Field<Type>&) {}

        virtual void operator*=(const Field<scalar>&) {}
        virtual void operator/=(const Field<scalar>&) {}

        virtual void operator=(const Type&) {}
        virtual void operator+=(const Type&) {}
        virtual void operator-=(const Type&) {}
        virtual void operator*=(const scalar) {}
        virtual void operator/=(const scalar) {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
