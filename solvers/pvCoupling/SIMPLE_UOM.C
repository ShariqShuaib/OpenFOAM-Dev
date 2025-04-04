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
    SIMPLE_UOM

Description
    We implement a simple implementation of the std SIMPLE algoritm present in 
    in "Versteeg, H.K. and Malalasekera, W. (2007) An Introduction to Computational Fluid Dynamics"..
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
        fvVectorMatrix UEqn
        (   
            // div(U * phi) - div (nu * grad(phi)) = -div(p) , Note in this case phi is u,v,w scalar components of velocities
            fvm::div(phi_star,U_star) - fvm::laplacian(nu,U_star) == -fvc::grad(p_star) 
        );
        // Solve for U*
        UEqn.solve();

        // a_p * U_p' = sum(i->n){a_i * U_i'} + sources - div(p')
        // U_p' = sum(i->n){a_i * U_i'} / a_p + sources / a_p - div(p') / a_p
        // For Std SIMPLE ignore contri from neighbours
        // i.e. U_p' = div(p') / a_p
        volScalarField A = UEqn.A();
        surfaceScalarField A_inv = fvc::interpolate(1/A);
        phi_star = fvc::interpolate(U_star) & mesh.Sf();
        

        // Now div(U) = 0
        // div(U* + U') = 0
        // div(U*) + div(U') = 0
        // div(U') = -div(U*)
        // div(-1/a_p * div(p')) = -div(U*)
        fvScalarMatrix pEqn
        (
            fvm::laplacian(A_inv,p_cor) == fvc::div(phi_star)
        );

        // Solve for p'
        pEqn.setReference(pRefCell,pRefValue);
        pEqn.solve();

        // Correct p as p = p* + p'*alpha with relaxation
        p = p_star + (alpha * p_cor);

        // Correct U as U = U* + U' , with U' = div(-1/a_p * div(p')) also apply relaxation
        U = U_star - 1/A * (fvc::grad(p_cor));
        U = alpha * U + (1 - alpha) * U_old;

        // Correct phi now
        phi = fvc::interpolate(U) & mesh.Sf();

        // Set p* = p ; U* = U ; phi* = phi and repeat
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
