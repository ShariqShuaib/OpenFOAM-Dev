/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) YEAR AUTHOR,AFFILIATION
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

Application
    SIMPLER

Description
    We implement a simple implementation of the std SIMPLER algoritm as presented
    in Malashk.
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"

    // Read-in the relaxation factor
    scalar alpha;
    fvSolution.lookup("alpha") >> alpha;


    // Read-in the Reference cell for pressure
    scalar pRefCell;
    fvSolution.lookup("pRefCell") >> pRefCell;

    // Read-in the Reference cell value
    scalar pRefValue;
    fvSolution.lookup("pRefValue") >> pRefValue;


    scalar cumulativeContErr = 0; 

    // Initilise U* , p* , phi* and p' NOTE: U = U* + U' and so on.
    volVectorField U_star(U);
    volScalarField p_star(p);
    surfaceScalarField phi_star(phi);
    volScalarField p_cor(p);

    while(runTime.loop())
    {
        Info << nl << "Iteration: " << runTime.timeName() << endl;
        #include "continuityErrs.H"

        // Note : WE DONT SOLVE THIS EQN JUST GET THE MATRIX OF COEFF TO SOLVE PRESSURE
        fvVectorMatrix UEqn
        (   
            // div(U * phi) - div (nu * grad(phi)) = -div(p) , Note in this case phi is u,v,w scalar components of velocities
            fvm::div(phi_star,U_star) - fvm::laplacian(nu,U_star) == -fvc::grad(p_star) 
        );

        volScalarField A = UEqn.A();
        surfaceScalarField A_inv = fvc::interpolate(1/A);
        volVectorField HbyA = UEqn.H() / A;

        // Solve the pressure equation, include the neighbouring corrections
        // Note this eqn is solved for exact value of Pressure i.e. p not p' or p*
        fvScalarMatrix pEqn
        (
            fvm::laplacian(A_inv,p) == fvc::div(HbyA) + fvc::div(phi_star)
        );
        pEqn.setReference(pRefCell,pRefValue);
        pEqn.solve();

        // Add Relaxation
        p = p_star + (alpha * p_cor);


        // Now solve the momentum eqn to get U*, we use exact pressure here
        fvVectorMatrix UEqn1
        (   
            // div(U * phi) - div (nu * grad(phi)) = -div(p)
            fvm::div(phi_star,U_star) - fvm::laplacian(nu,U_star) == -fvc::grad(p) 
        );

        // Solve for U*
        UEqn1.solve();
        // Compute phi*
        phi_star = fvc::interpolate(U_star) & mesh.Sf();
        surfaceScalarField A_inv1 = fvc::interpolate(1/UEqn1.A());

        // We then solve pressure correction equation to compute the p' that must be added to 
        // U* to correct it
        fvScalarMatrix pEqn1
        (
            fvm::laplacian(A_inv1,p_cor) == fvc::div(phi_star)
        );
        pEqn1.setReference(pRefCell,pRefValue);
        
        // Solve for p'
        pEqn1.solve();

        // Correct U as U = U* + U' , with U' = div(-1/a_p * div(p')) also apply relaxation
        U = U_star - 1/A * (fvc::grad(p_cor));
        U = alpha * U + (1 - alpha) * U_old;

        // Correct phi now
        phi = fvc::interpolate(U) & mesh.Sf();

        // Set p* = p ; U* = U ; phi* = phi and repeat
        // Note in this case we ideally dont have p_star its same as p so you can ignore this or replace this
        p_star = p;
        U_star = U;
        phi_star = phi;


        // Correct BCs
        U.correctBoundaryConditions();
        p.correctBoundaryConditions();

        // Save for relaxation
        U_old = U;

        runTime.write();
    }
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< nl;
    runTime.printExecutionTime(Info);

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
