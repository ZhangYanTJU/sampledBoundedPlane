/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
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

#include "sampledBoundedPlane.H"
#include "dictionary.H"
#include "polyMesh.H"
#include "volFields.H"

#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(sampledBoundedPlane, 0);
    addNamedToRunTimeSelectionTable(sampledSurface, sampledBoundedPlane, word, boundedPlane);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sampledBoundedPlane::sampledBoundedPlane
(
    const word& name,
    const polyMesh& mesh,
    const plane& planeDesc,
    const keyType& zoneKey,
    const bool triangulate
)
:
    sampledSurface(name, mesh),
    cuttingPlane(planeDesc),
    zoneKey_(zoneKey),
    bounds_(),
    triangulate_(triangulate),
    needsUpdate_(true)
{
    if (debug && zoneKey_.size() && mesh.cellZones().findIndex(zoneKey_) < 0)
    {
        Info<< "cellZone " << zoneKey_
            << " not found - using entire mesh" << endl;
    }
}


Foam::sampledBoundedPlane::sampledBoundedPlane
(
    const word& name,
    const polyMesh& mesh,
    const dictionary& dict
)
:
    sampledSurface(name, mesh, dict),
    cuttingPlane(plane(dict)),
    zoneKey_(keyType::null),
    bounds_(dict.lookupOrDefault("bounds", boundBox::invertedBox)),
    triangulate_(dict.lookupOrDefault("triangulate", true)),
    needsUpdate_(true)
{
    // Make plane relative to the coordinateSystem (Cartesian)
    // allow lookup from global coordinate systems
    if (dict.found("coordinateSystem"))
    {
        coordinateSystem cs(mesh, dict.subDict("coordinateSystem"));

        point  base = cs.globalPosition(planeDesc().refPoint());
        vector norm = cs.globalVector(planeDesc().normal());

        // Assign the plane description
        static_cast<plane&>(*this) = plane(base, norm);
    }

    dict.readIfPresent("zone", zoneKey_);

    if (debug && zoneKey_.size() && mesh.cellZones().findIndex(zoneKey_) < 0)
    {
        Info<< "cellZone " << zoneKey_
            << " not found - using entire mesh" << endl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::sampledBoundedPlane::~sampledBoundedPlane()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::sampledBoundedPlane::needsUpdate() const
{
    return needsUpdate_;
}


bool Foam::sampledBoundedPlane::expire()
{
    // Already marked as expired
    if (needsUpdate_)
    {
        return false;
    }

    sampledSurface::clearGeom();

    needsUpdate_ = true;
    return true;
}


bool Foam::sampledBoundedPlane::update()
{
    Info<<"I am here =========sampledBoundedPlane======="<<endl;
    if (!needsUpdate_)
    {
        return false;
    }

    sampledSurface::clearGeom();

    labelList selectedCells = mesh().cellZones().findMatching(zoneKey_).used();
    
    bool fullMesh = returnReduce(selectedCells.empty(), andOp<bool>());

    //if (!bounds_.empty())
    //{
        const auto& cellCentres = static_cast<const fvMesh&>(mesh()).C();

        if (fullMesh)
        {
            const label len = mesh().nCells();

            selectedCells.setSize(len);

            label count = 0;
            for (label celli=0; celli < len; ++celli)
            {
                if (bounds_.contains(cellCentres[celli]))
                {
                    selectedCells[count++] = celli;
                }
            }

            selectedCells.setSize(count);
        }
        else
        {
            label count = 0;
            for (const label celli : selectedCells)
            {
                if (bounds_.contains(cellCentres[celli]))
                {
                    selectedCells[count++] = celli;
                }
            }

            selectedCells.setSize(count);
        }

        fullMesh = false;
    //}
    
    //if (selectedCells.empty())
    if (fullMesh)
    {
        reCut(mesh(), triangulate_);
    }
    else
    {
        reCut(mesh(), triangulate_, selectedCells);
    }

    if (debug)
    {
        print(Pout);
        Pout<< endl;
    }

    needsUpdate_ = false;
    return true;
}


Foam::tmp<Foam::scalarField> Foam::sampledBoundedPlane::sample
(
    const volScalarField& vField
) const
{
    return sampleField(vField);
}


Foam::tmp<Foam::vectorField> Foam::sampledBoundedPlane::sample
(
    const volVectorField& vField
) const
{
    return sampleField(vField);
}


Foam::tmp<Foam::sphericalTensorField> Foam::sampledBoundedPlane::sample
(
    const volSphericalTensorField& vField
) const
{
    return sampleField(vField);
}


Foam::tmp<Foam::symmTensorField> Foam::sampledBoundedPlane::sample
(
    const volSymmTensorField& vField
) const
{
    return sampleField(vField);
}


Foam::tmp<Foam::tensorField> Foam::sampledBoundedPlane::sample
(
    const volTensorField& vField
) const
{
    return sampleField(vField);
}


Foam::tmp<Foam::scalarField> Foam::sampledBoundedPlane::interpolate
(
    const interpolation<scalar>& interpolator
) const
{
    return interpolateField(interpolator);
}


Foam::tmp<Foam::vectorField> Foam::sampledBoundedPlane::interpolate
(
    const interpolation<vector>& interpolator
) const
{
    return interpolateField(interpolator);
}

Foam::tmp<Foam::sphericalTensorField> Foam::sampledBoundedPlane::interpolate
(
    const interpolation<sphericalTensor>& interpolator
) const
{
    return interpolateField(interpolator);
}


Foam::tmp<Foam::symmTensorField> Foam::sampledBoundedPlane::interpolate
(
    const interpolation<symmTensor>& interpolator
) const
{
    return interpolateField(interpolator);
}


Foam::tmp<Foam::tensorField> Foam::sampledBoundedPlane::interpolate
(
    const interpolation<tensor>& interpolator
) const
{
    return interpolateField(interpolator);
}


void Foam::sampledBoundedPlane::print(Ostream& os) const
{
    os  << "sampledBoundedPlane: " << name() << " :"
        << "  base:" << refPoint()
        << "  normal:" << normal()
        << "  triangulate:" << triangulate_
        << "  faces:" << faces().size()
        << "  points:" << points().size();
}


// ************************************************************************* //
