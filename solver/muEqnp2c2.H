
//! The unit normal vector to the interface with a small number in denominator to prevent solution from diverging
volVectorField n=dimx*fvc::grad(phi_1)/(1E-20+mag(dimx*fvc::grad(phi_1)));

//! Chemical potential equation with approximate slope of c-mu curve according to the parabolic approximation for free energy
if (phases == 2)
{
hphi = phi_1*phi_1*phi_1*(6*phi_1*phi_1 - 15*phi_1 +10);
if (components == 2)
{
fvScalarMatrix muEqn
(
    (0.5)*((hphi/A_Sol) + ((1-hphi)/A_Liq))*dimt*fvm::ddt(mu_1) == dimx*dimx*fvm::laplacian(diff_Sol*0.5*(1/A_Sol)*phi_1 + diff_Liq*0.5*(1/A_Liq)*(1-phi_1), mu_1) - 0.5*((mu_1 - B_Sol)/A_Sol - (mu_1 - B_Liq)/A_Liq)*dimt*fvc::ddt(phi_1)*30.0*phi_1*phi_1*(1.0-phi_1)*(1.0-phi_1) - Tdot*(((-dBdT)/(2*A_Sol) + (mu_1 - B_Sol)*A_SoldT/(-2*A_Sol*A_Sol))*hphi + (mu_1 - B_Liq)*A_LiqdT*(1-hphi)/(-2*A_Liq*A_Liq)) - anti_trap*epsilon*0.5*((mu_1 - B_Sol)/A_Sol - (mu_1 - B_Liq)/A_Liq)*dimx*fvc::div((n*dimt*fvc::ddt(phi_1)))
);

//muEqn.solve();

InitialResidual_10 = muEqn.solve().max().initialResidual();

    Info<< "Min/max mu:" << min(mu_1).value() << ' ' << max(mu_1).value() << endl;
}

}

