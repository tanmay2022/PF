
//! The unit normal vector to the interface with a small number in denominator to prevent solution from diverging
volVectorField n=dimx*fvc::grad(phi_1)/(1E-20+mag(dimx*fvc::grad(phi_1)));

volScalarField hphi1 = 0; //-phi*phi*phi*(10*phi*phi*phi - 36*phi*phi + 45*phi - 20);
volScalarField hphi2 = 0;
volScalarField hphi3 = 0;
volScalarField hphi4 = 0;

//! Chemical potential equation with approximate slope of c-mu curve according to the parabolic approximation for free energy
if (phases == 4)
{
hphi1 = phi_1*phi_1*(3-2*phi_1) + 2*phi_1*(phi_2*phi_3 + phi_2*phi_4 + phi_3*phi_4);
hphi2 = phi_2*phi_2*(3-2*phi_2) + 2*phi_2*(phi_1*phi_3 + phi_1*phi_4 + phi_3*phi_4);
hphi3 = phi_3*phi_3*(3-2*phi_3) + 2*phi_3*(phi_1*phi_2 + phi_1*phi_4 + phi_2*phi_4);
hphi4 = phi_4*phi_4*(3-2*phi_4) + 2*phi_4*(phi_1*phi_2 + phi_1*phi_3 + phi_2*phi_3);
if (components == 2)
{
fvScalarMatrix muEqn
(
    (0.5)*((hphi1/A_Sol1) + (hphi2/A_Sol2) + (hphi3/A_Sol3) + (hphi4/A_Liq))*dimt*fvm::ddt(mu_1) == dimx*dimx*fvm::laplacian(diff_Sol*0.5*(1/A_Sol)*phi_1 + diff_Liq*0.5*(1/A_Liq)*(1-phi_1), mu_1) 
    - 0.5*((mu_1 - B_Sol1)*dimt*fvc::ddt(phi_1)*(6.0*phi_1*(1-phi_1) + 2*(phi_2*phi_3 + phi_2*phi_4 + phi_3*phi_4))/A_Sol1 + (mu_1 - B_Sol2)*dimt*fvc::ddt(phi_2)*(6.0*phi_2*(1-phi_2) + 2*(phi_1*phi_3 + phi_1*phi_4 + phi_3*phi_4))/A_Sol2 + (mu_1 - B_Sol3)*dimt*fvc::ddt(phi_3)*(6.0*phi_3*(1-phi_3) + 2*(phi_1*phi_2 + phi_1*phi_4 + phi_2*phi_4))/A_Sol3 + (mu_1 - B_Liq)*dimt*fvc::ddt(phi_4)*(6.0*phi_4*(1-phi_4) + 2*(phi_1*phi_2 + phi_1*phi_3 + phi_2*phi_3))/A_Liq) 
    - Tdot*(((-dB1dT)/(2*A_Sol1) + (mu_1 - B_Sol1)*A_Sol1dT/(-2*A_Sol1*A_Sol1))*hphi1 + ((-dB2dT)/(2*A_Sol2) + (mu_1 - B_Sol2)*A_Sol2dT/(-2*A_Sol2*A_Sol2))*hphi2 + ((-dB3dT)/(2*A_Sol3) + (mu_1 - B_Sol3)*A_Sol3dT/(-2*A_Sol3*A_Sol3))*hphi3 + (mu_1 - B_Liq)*A_LiqdT*hphi4/(-2*A_Liq*A_Liq))
     //- anti_trap*epsilon*0.5*((mu_1 - B_Sol)/A_Sol - (mu_1 - B_Liq)/A_Liq)*dimx*fvc::div((n*dimt*fvc::ddt(phi_1)))
);

//muEqn.solve();

InitialResidual_10 = muEqn.solve().max().initialResidual();

    Info<< "Min/max mu:" << min(mu_1).value() << ' ' << max(mu_1).value() << endl;
}

}

