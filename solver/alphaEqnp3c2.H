#include "dAdgradPhiMod1.H"


if (phases == 3)
{
if (components == 2)
{
//For phi_1
    
    fvScalarMatrix phi_1Eqn
    (
        dimt*omega*epsilon*fvm::ddt(phi_1) -2.0*gamma_1*epsilon*dimx*dimx*fvm::laplacian(DTs, phi_1) == -(18.0*gamma_1*phi_1*(1.0 -phi_1)*(1.0 -2.0*phi_1))/(epsilon)
        - ((-1.0/4.0*(mu_1 - B_Sol)*(mu_1 - B_Sol)/A_Sol + D_Sol)*0.25*15*(3.0*phi_1*phi_1*phi_1*phi_1 - 4.0*phi_1*phi_1*phi_1 - phi_1*phi_1 + 2.0*phi_1 - (phi_2 - phi_3)*(phi_2 - phi_3)*(2.0*phi_1 - 3.0*phi_1*phi_1))
         + (-1.0/4.0*(mu_1 - B_Sol)*(mu_1 - B_Sol)/A_Sol + D_Sol)*-0.5*15*((phi_2*phi_2 - phi_2*phi_2*phi_2)*(phi_1 - phi_3))
          + (-1.0/4.0*(mu_1 - B_Liq)*(mu_1 - B_Liq)/A_Liq + D_Liq)*-0.5*15*((phi_3*phi_3 - phi_3*phi_3*phi_3)*(phi_1 - phi_2)))*(1 + noise_mag*randNum.globalScalar01())/Vm
             - (6.0*phi_1*(1-phi_1)*(0.5*(deltaSigmaD && (symm(fvc::grad(D))-phi_1*phi_1*(3-2*phi_1)*cEigenStrain1 - phi_2*phi_2*(3-2*phi_2)*cEigenStrain2
                                                                    ))-(Sigma && cEigenStrain1)))
    );
//    phi_1Eqn -= 6.0*phi_1*(1-phi_1)*(0.5*(deltaSigmaD && (symm(fvc::grad(D))-phi_1*phi_1*(3-2*phi_1)*cEigenStrain1 - phi_2*phi_2*(3-2*phi_2)*cEigenStrain2 
//								    - phi_3*phi_3*(3-2*phi_3)*cEigenStrain3))-(Sigma && cEigenStrain1));

    
//For phi_2        
    
    fvScalarMatrix phi_2Eqn
    (
        dimt*omega*epsilon*fvm::ddt(phi_2) -2.0*gamma_2*epsilon*dimx*dimx*fvm::laplacian(DTs, phi_2) == -(18.0*gamma_2*phi_2*(1.0 -phi_2)*(1.0 -2.0*phi_2))/(epsilon)
        - ((-1.0/4.0*(mu_1 - B_Sol)*(mu_1 - B_Sol)/A_Sol + D_Sol)*-0.5*15*((phi_1*phi_1 - phi_1*phi_1*phi_1)*(phi_2 - phi_3))
         + (-1.0/4.0*(mu_1 - B_Sol)*(mu_1 - B_Sol)/A_Sol + D_Sol)*0.25*15*(3.0*phi_2*phi_2*phi_2*phi_2 - 4.0*phi_2*phi_2*phi_2 - phi_2*phi_2 + 2.0*phi_2 - (phi_3 - phi_1)*(phi_3 - phi_1)*(2.0*phi_2 - 3.0*phi_2*phi_2))
          + (-1.0/4.0*(mu_1 - B_Liq)*(mu_1 - B_Liq)/A_Liq + D_Liq)*-0.5*15*((phi_3*phi_3 - phi_3*phi_3*phi_3)*(phi_2 - phi_1)))*(1 + noise_mag*randNum.globalScalar01())/Vm
		- (6.0*phi_2*(1-phi_2)*(0.5*(deltaSigmaD && (symm(fvc::grad(D))-phi_1*phi_1*(3-2*phi_1)*cEigenStrain1 - phi_2*phi_2*(3-2*phi_2)*cEigenStrain2
                                                           ))-(Sigma && cEigenStrain2)))
    );
//    phi_2Eqn -= 6.0*phi_2*(1-phi_2)*(0.5*(deltaSigmaD && (symm(fvc::grad(D))-phi_1*phi_1*(3-2*phi_1)*cEigenStrain1 - phi_2*phi_2*(3-2*phi_2)*cEigenStrain2 
//								    - phi_3*phi_3*(3-2*phi_3)*cEigenStrain3))-(Sigma && cEigenStrain2));

    
//For phi_3  

     deltaSigmaD = ((mu2_elast)*twoSymm(fvc::grad(D)) + (lambda2)*(I*tr(fvc::grad(D))) 
		   - 2*(mu2_elast)*(phi_1*phi_1*(3-2*phi_1)*cEigenStrain1 + phi_2*phi_2*(3-2*phi_2)*cEigenStrain2) 
		   - (lambda2)*(I*tr((phi_1*phi_1*(3-2*phi_1)*cEigenStrain1 + phi_2*phi_2*(3-2*phi_2)*cEigenStrain2))))
		   + (mu2_elast_)*strain;
    
    fvScalarMatrix phi_3Eqn
    (
        dimt*omega*epsilon*fvm::ddt(phi_3) -2.0*gamma_3*epsilon*dimx*dimx*fvm::laplacian(DTs, phi_3) == -(18.0*gamma_3*phi_3*(1.0 -phi_3)*(1.0 -2.0*phi_3))/(epsilon)
        - ((-1.0/4.0*(mu_1 - B_Sol)*(mu_1 - B_Sol)/A_Sol + D_Sol)*-0.5*15*((phi_1*phi_1 - phi_1*phi_1*phi_1)*(phi_3 - phi_2))
         + (-1.0/4.0*(mu_1 - B_Sol)*(mu_1 - B_Sol)/A_Sol + D_Sol)*-0.5*15*((phi_2*phi_2 - phi_2*phi_2*phi_2)*(phi_3 - phi_1))
          + (-1.0/4.0*(mu_1 - B_Liq)*(mu_1 - B_Liq)/A_Liq + D_Liq)*0.25*15*(3.0*phi_3*phi_3*phi_3*phi_3 - 4.0*phi_3*phi_3*phi_3 - phi_3*phi_3 + 2.0*phi_3 - (phi_1 - phi_2)*(phi_1 - phi_2)*(2.0*phi_3 - 3.0*phi_3*phi_3)))*(1 + noise_mag*randNum.globalScalar01())/Vm
			- (6.0*phi_3*(1-phi_3)*(0.5*(deltaSigmaD && (symm(fvc::grad(D))-phi_1*phi_1*(3-2*phi_1)*cEigenStrain1 - phi_2*phi_2*(3-2*phi_2)*cEigenStrain2))))
    );
//    phi_3Eqn -= 6.0*phi_3*(1-phi_3)*(0.5*(deltaSigmaD && (symm(fvc::grad(D))-phi_1*phi_1*(3-2*phi_1)*cEigenStrain1 - phi_2*phi_2*(3-2*phi_2)*cEigenStrain2 
//								    - phi_3*phi_3*(3-2*phi_3)*cEigenStrain3))-(Sigma && cEigenStrain3));    

    
    
    
    phi_1Eqn.solve();
    phi_2Eqn.solve();
    phi_3Eqn.solve();

  Info<< "Min/max sigmaproduct:" << min((Sigma && cEigenStrain2)()).value() << '\t'<< max((Sigma && cEigenStrain2)()).value() << endl;
 
//    scalar numphases = 3.0;   
    
    volScalarField lambda_phase = dimt*((fvc::ddt(phi_1) + fvc::ddt(phi_2) + fvc::ddt(phi_3))/phases);
    Info<< "Min/max lambda_phase:" << min(lambda_phase()).value() << '\t'<< max(lambda_phase()).value() << endl;


//=============== Substracting lambda_phase from evolution equations ==============================//  
    
  //Update the phase-values such that sum upto 1
    phi_1 -=  lambda_phase*runTime.deltaTValue();
    phi_2 -=  lambda_phase*runTime.deltaTValue();
    phi_3 -=  lambda_phase*runTime.deltaTValue();
}
}

