/*---------------------------------------------------------------------------*\
Copyright (C) 2022 Matthias Renaud, Cyrille Bonamy, Julien Chauchat
                   and contributors

License
    This file is part of ScourFOAM.

    ScourFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    ScourFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with ScourFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/
volScalarField rAU(1.0/UEqn.A());
surfaceScalarField rAUf("rAUf", fvc::interpolate(rAU));
volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p_rgh));

surfaceScalarField phig(-rAUf*ghf*fvc::snGrad(rhok)*mesh.magSf());

surfaceScalarField phiHbyA("phiHbyA", fvc::flux(HbyA));

if (pimple.ddtCorr())
{
    phiHbyA += MRF.zeroFilter(rAUf*fvc::ddtCorr(U, phi, Uf));
}

// Add Boussinesq term 
//if (switchSuspension=="on")
//{
  phiHbyA += phig;
//}

MRF.makeRelative(phiHbyA);

if (p.needReference())
{
    fvc::makeRelative(phiHbyA, U);
    adjustPhi(phiHbyA, U, p_rgh);
    fvc::makeAbsolute(phiHbyA, U);
}

tmp<volScalarField> rAtU(rAU);

if (pimple.consistent())
{
    rAtU = 1.0/max(1.0/rAU - UEqn.H1(), 0.1/rAU);
    phiHbyA +=
        fvc::interpolate(rAtU() - rAU)*fvc::snGrad(p_rgh)*mesh.magSf();
    HbyA -= (rAU - rAtU())*fvc::grad(p_rgh);
}

if (pimple.nCorrPISO() <= 1)
{
    tUEqn.clear();
}

// Update the pressure BCs to ensure flux consistency
constrainPressure(p, U, phiHbyA, rAtU(), MRF);

// Non-orthogonal pressure corrector loop
while (pimple.correctNonOrthogonal())
{
    fvScalarMatrix p_rghEqn
    (
        fvm::laplacian(rAtU(), p_rgh) == fvc::div(phiHbyA)
    );

    p_rghEqn.setReference(pRefCell, pRefValue);

    p_rghEqn.solve(mesh.solver(p_rgh.select(pimple.finalInnerIter())));

    if (pimple.finalNonOrthogonalIter())
    {
        phi = phiHbyA - p_rghEqn.flux();

        // Explicitly relax pressure for momentum corrector
        p_rgh.relax();

        // Add Boussinesq term 
        //if (switchSuspension=="on")
        //{
        U = HbyA + rAU*fvc::reconstruct((phig - p_rghEqn.flux())/rAUf);
        //}
        //else
        // {
        //    U = HbyA - rAtU*fvc::grad(p_rgh);
        //}   
        
        U.correctBoundaryConditions();
        fvOptions.correct(U);
        
        // Correct Uf if the mesh is moving
        fvc::correctUf(Uf, U, phi);

        // Make the fluxes relative to the mesh motion
        fvc::makeRelative(phi, U);
    }
}

#include "continuityErrs.H"

p = p_rgh + rhok*gh;

if (p_rgh.needReference())
{
    p += dimensionedScalar
        (
            "p",
            p.dimensions(),
            pRefValue - getRefCellValue(p, pRefCell)
        );
    p_rgh = p - rhok*gh;
}
/* 
// Explicitly relax pressure for momentum corrector
p.relax();

U = HbyA - rAtU*fvc::grad(p);
U.correctBoundaryConditions();
fvOptions.correct(U);

//phip = fvc::interpolate(U+Ws) & mesh.Sf();

// Correct Uf if the mesh is moving
fvc::correctUf(Uf, U, phi);

// Make the fluxes relative to the mesh motion
fvc::makeRelative(phi, U);
*/
