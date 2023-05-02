
//! The unit normal vector to the interface with a small number in denominator to prevent solution from diverging
volVectorField n=dimx*fvc::grad(phi)/(1E-20+mag(dimx*fvc::grad(phi)));

volScalarField hphi = phi*phi*phi*(6*phi*phi - 15*phi +10); //-phi*phi*phi*(10*phi*phi*phi - 36*phi*phi + 45*phi - 20);

//! Chemical potential equation with approximate slope of c-mu curve according to the parabolic approximation for free energy
fvScalarMatrix muEqn
(
    (0.5)*((hphi/A_Sol) + ((1-hphi)/A_Liq))*dimt*fvm::ddt(mu) == dimx*dimx*fvm::laplacian(diff_Sol*0.5*(1/A_Sol)*phi + diff_Liq*0.5*(1/A_Liq)*(1-phi), mu) - 0.5*((mu - B_Sol)/A_Sol - (mu - B_Liq)/A_Liq)*dimt*fvc::ddt(phi)*30.0*phi*phi*(1.0-phi)*(1.0-phi) - Tdot*(((-dBdT)/(2*A_Sol) + (mu - B_Sol)*A_SoldT/(-2*A_Sol*A_Sol))*hphi + (mu - B_Liq)*A_LiqdT*(1-hphi)/(-2*A_Liq*A_Liq)) - anti_trap*epsilon*0.5*((mu - B_Sol)/A_Sol - (mu - B_Liq)/A_Liq)*dimx*fvc::div((n*dimt*fvc::ddt(phi)))
);

//muEqn.solve();

InitialResidual_10 = muEqn.solve().max().initialResidual();

    Info<< "Min/max mu:" << min(mu).value() << ' ' << max(mu).value() << endl;

