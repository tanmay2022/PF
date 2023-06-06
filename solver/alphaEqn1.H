

//! Implicit discretization using fvm class for time derivative and laplacian.
//! Explicit discretization using fvc class for divergence.
//! Phase-field equation with approximate relation between c, mu and T according to the parabolic approximation for free energy
if (phases == 2)
{

if (components == 3)
{
fvScalarMatrix phi_1Eqn (
/*omega*epsilon*dimt*fvm::ddt(alpha)
	-epsilon*gamma*dimx*dimx*fvm::laplacian(alpha)
	+2.0*alpha*(alpha-1.0)*(2.0*alpha-1.0)/epsilon		*/

omega*epsilon*dimt*fvm::ddt(phi_1)     //alpha=1 liquid phase
                ==
// 		2.0*epsilon*gamma*dimx*dimx*fvm::laplacian(ac_01*ac_01,phi) + 2*gamma*epsilon*dimx*fvc::div(dadgradPhi) - 18.0*(gamma/epsilon)*(phi)*(1-phi)*(1-2.0*phi) 
//                 + ((ceq_1[0][1] - ceq_2[0][1])*(mu_1 - mu_eq[0][1]) + (ceq_1[1][1] - ceq_2[1][1])*(mu_2 - mu_eq[1][1]))*30.0*phi*phi*(1.0-phi)*(1.0-phi)
    2.0*epsilon*gamma_1*dimx*dimx*fvm::laplacian(ac_01*ac_01,phi_1) + 2*gamma_1*epsilon*dimx*fvc::div(dadgradPhi) - 18.0*(gamma_1/epsilon)*(phi_1)*(1-phi_1)*(1-2.0*phi_1) 
    + ((-0.5*dmudc_l[0][0]*(dcdmu_l[0][0]*mu_1 + dcdmu_l[0][1]*mu_2)*(dcdmu_l[0][0]*mu_1 + dcdmu_l[0][1]*mu_2) - dmudc_l[0][1]*(dcdmu_l[0][0]*mu_1 + dcdmu_l[0][1]*mu_2)*(dcdmu_l[1][0]*mu_1 + dcdmu_l[1][1]*mu_2) - 0.5*dmudc_l[1][1]*(dcdmu_l[1][0]*mu_1 + dcdmu_l[1][1]*mu_2)*(dcdmu_l[1][0]*mu_1 + dcdmu_l[1][1]*mu_2)) - (-0.5*dmudc_a[0][0]*(dcdmu_a[0][0]*(mu_1 - B_a1) + dcdmu_a[0][1]*(mu_2 - B_a2))*(dcdmu_a[0][0]*(mu_1 - B_a1) + dcdmu_a[0][1]*(mu_2 - B_a2)) - dmudc_a[0][1]*(dcdmu_a[0][0]*(mu_1 - B_a1) + dcdmu_a[0][1]*(mu_2 - B_a2))*(dcdmu_a[1][0]*(mu_1 - B_a1) + dcdmu_a[1][1]*(mu_2 - B_a2)) - 0.5*dmudc_a[1][1]*(dcdmu_a[1][0]*(mu_1 - B_a1) + dcdmu_a[1][1]*(mu_2 - B_a2))*(dcdmu_a[1][0]*(mu_1 - B_a1) + dcdmu_a[1][1]*(mu_2 - B_a2)) + DD_a))*30.0*phi_1*phi_1*(1.0-phi_1)*(1.0-phi_1)*(1 + noise_mag*randNum.globalScalar01())/Vm - epsilon*mu_th*mu_th*grad_qt_sqr*((3*(7*phi_1*phi_1+8*phi_1*phi_1*phi_1)*(1-phi_1) + 2*phi_1*phi_1*phi_1*(7-6*phi_1))/(1+eps-phi_1))
    - 6.0*phi_1*(1-phi_1)*(0.5*(deltaSigmaD && (symm(fvc::grad(D))-phi_1*phi_1*(3-2*phi_1)*cEigenStrain1))-(Sigma && cEigenStrain1))
);

//alphaEqn.solve();

InitialResidual_0 = phi_1Eqn.solve().max().initialResidual();
}
}


