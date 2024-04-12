/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2024 2024 Konstantinos Missios, Roskilde University
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

#include "interfaceArea.H"
#include "Time.H"
#include "fvMesh.H"
#include "fvcGrad.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(interfaceArea, 0);
    addToRunTimeSelectionTable(functionObject, interfaceArea, dictionary);
}
}

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::interfaceArea::writeFileHeader(Ostream& os)
{
    writeHeader(os, "Interface area evolution");
    writeCommented(os, "Time");
    writeCommented(os, "Interface area");

    os  << endl;
}



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::interfaceArea::interfaceArea
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    writeFile(mesh_, name, typeName, dict),
    // boolData_(dict.getOrDefault<bool>("boolData"), true),
    // labelData_(dict.get<label>("labelData")),
    phaseName_(dict.getOrDefault<word>("phaseName", "phase1"))
    // scalarData_(dict.getOrDefault<scalar>("scalarData", 1.0))
{
    // read(dict);    
    if (read(dict))
    {
        writeFileHeader(file());
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::interfaceArea::read(const dictionary& dict)
{
    dict.readIfPresent("phaseName", phaseName_);

    return true;
}


bool Foam::functionObjects::interfaceArea::execute()
{
    volScalarField alpha =  mesh_.lookupObject<volScalarField>("alpha." + phaseName_);
            
    volVectorField gradAlpha(fvc::grad(alpha));       
     
    volScalarField magGradAlpha = mag(gradAlpha);

    interfaceArea_ = gSum(magGradAlpha*mesh_.V().field());
    
    return true;
}



bool Foam::functionObjects::interfaceArea::write()
{
    Ostream& os = file();
    
    writeCurrentTime(os);
    
    os  << interfaceArea_ << endl;

    return true;
}


// ************************************************************************* //
