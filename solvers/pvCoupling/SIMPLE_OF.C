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
    SIMPLE_OF

Description
    We implement a simple implementation of SIMPLE algoritm this is OF implem
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

    while(runTime.loop())
    {
        Info << nl << "Iteration: " << runTime.timeName() << endl;

        fvVectorMatrix UEqn
        (
            // div(U * phi) - div (nu * grad(phi)) = -div(p) , Note in this case phi is u,v,w scalar components of velocities
            fvm::div(phi,U) - fvm::laplacian(nu,U) == -fvc::grad(p)
        );

        // Solving for U with initi values of pressure , note these will not satify the continity
        UEqn.solve(); 
    
        // a_p * U_p = sum(i->n){a_i * U_i} + sources - div(p) , where n total no of adj neighbours
        // Now U_p = sum(i->n){a_i * U_i} / a_p - div(p) / a_p , note sources are taken in um(i->n){a_i * U_i}
        // eg [1 2 1 0 0][U1 U_p U3 U4 U5] = [Sink] -> rearrange for U_p 
        // 2 * U_p = -U1 - U3 + Sink + pressure 
        // Now here 2 is A matrix and -U1 - U3 + **Sink is the H ; H is vector as it has 2 compons if solving u,v velo eqns

        volScalarField A = UEqn.A();  // This is U_p a diagonal matrix
        volVectorField H = UEqn.H();  // This is a_i * U_i matrix
        volVectorField HbyA = H / A;
        

        // div(U) = 0
        // div(H / A , -1/A * div(P)) = 0
        // solving 
        // div(1/A * div(p)) = div(H/A)
        // solve for p
        // to use fvm :: 1/A has to be flux i.e. fvm
        
        surfaceScalarField A_inv = fvc::interpolate(1/A);
        fvScalarMatrix pEqn
        (
            fvm::laplacian(A_inv,p) == fvc::div(HbyA)
        );

        // Ref pressure i.e. abs
        pEqn.setReference(pRefCell,pRefValue);
        pEqn.solve();

        // Apply the tolerance
        p = alpha * p + (1 - alpha) * p_old;

        // Update the velocity
        U = HbyA - 1/A * (fvc::grad(p));

        phi = fvc::interpolate(U) & mesh.Sf();


        // Update the boundary conditions
        U.correctBoundaryConditions();
        p.correctBoundaryConditions();

        p_old = p;

        #include "continuityErrs.H"

        runTime.write();
    }
   

    


    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< nl;
    runTime.printExecutionTime(Info);

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
