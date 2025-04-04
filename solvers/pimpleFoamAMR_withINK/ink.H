singlePhaseTransportModel lamTrans(U, phi);
const volScalarField& nu = lamTrans.nu();

fvScalarMatrix inkEqn
(
   fvm::ddt(ink) + fvm::div(phi,ink) - fvm::laplacian(nu,ink)
);

inkEqn.solve();
ink.relax();