/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2024 Konstantinos Missios, Roskilde University
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
    Foam::functionObjects::interfaceArea

Group
    grpFieldFunctionObjects

Description
    Evaluated the interface area between two fluids based on the information 
    provided by the volume fraction districution. The evaluation is given as
    
    \f[
        A_i =\sum_{C=1}^N \lvert\nabla \alpha_C\rvert V_C
    \f]

    where
    \vartable
      A_i         | total interface area [m^2]
      C           | cell index [#]
      \alpha_C    | volume fraction of cell C [#]
      V_C         | volume of cell C [m^3]
    \endvartable

    Operands:
    \table
      Operand      | Type | Location
      input        |  -   | -
      output file  | dat  | $FOAM_CASE/postProcessing/\<FO\>/\<time\>/\<file\>
    \endtable

Usage
    Minimal example by using \c system/controlDict.functions:
    \verbatim
    interfaceArea1
    {
        // Mandatory entries (unmodifiable)
        type           interfaceArea;
        libs           (interfaceAreaFunctionObject);

        // Mandatory entries (runtime modifiable)
        ...

        // Mandatory (inherited) entries (unmodifiable)
        ...

        // Mandatory (inherited) entries (runtime unmodifiable)
        ...

        // Optional entries (unmodifiable)
        ...

        // Optional entries (runtime modifiable)
        boolData       <bool>;
        labelData      <label>;
        wordData       <word>;
        scalarData     <scalar>;

        // Optional (inherited) entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property   | Description                      | Type | Req'd | Dflt
      type       | Type name: interfaceArea        | word |  yes  | -
      libs       | Library name: interfaceAreaFunctionObject  <!--
             --> | word  | yes | -
      boolData   | <explanation>                    | bool |  yes  | -
      labelData  | <explanation>                    | label |  yes | -
      wordData   | <explanation>                    | word |  yes  | -
      scalarData | <explanation>                    | scalar |  no | 1.0
      wordListData  | <explanation>                 | wordList | yes | -
    \endtable

    Options for the \c ENTRY entry:
    \verbatim
        <option1>
        <option2> | <explanation>
      ...
    \endverbatim

    The inherited entries are elaborated in:
     - \link functionObject.H \endlink
     - \link fieldExpression.H \endlink
     - \link fieldsExpression.H \endlink
     - \link writeFile.H \endlink
     ...

    <if \c postProcess is applicable>
    Minimal example by using the \c postProcess utility:
    \verbatim
        postProcess -func interfaceArea
    \endverbatim

    <if \c postProcess is not applicable>
    Usage by the \c postProcess utility is not available.

Note
    - <note1>
    - <note2>
    ...

See also
    - Foam::functionObject
    - Foam::functionObjects::fvMeshFunctionObject
    - ExtendedCodeGuide::functionObjects::field::interfaceArea
    ...

SourceFiles
    interfaceArea.C
    interfaceAreaTEMPLATES.C
    ...

\*---------------------------------------------------------------------------*/

#ifndef interfaceArea_H
#define interfaceArea_H

#include "fvMeshFunctionObject.H"
#include "writeFile.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                   Class interfaceArea Declaration
\*---------------------------------------------------------------------------*/

class interfaceArea
:
    public fvMeshFunctionObject,
    public writeFile
{
protected:

        //- name of the phase
        word phaseName_;

        //- interface Area
        scalar interfaceArea_;

        //- Output file header information
        virtual void writeFileHeader(Ostream& os);

public:

    //- Runtime type information
    TypeName("interfaceArea");


    // Constructors

        //- Construct from Time and dictionary
        interfaceArea
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- No copy construct
        interfaceArea(const interfaceArea&) = delete;

        //- No copy assignment
        void operator=(const interfaceArea&) = delete;


    //- Destructor
    virtual ~interfaceArea() = default;


    // Member Functions

        //- Read the interfaceArea data
        virtual bool read(const dictionary& dict);

        //- Calculate the interface are
        virtual bool execute();

        //- Write the interface area
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
