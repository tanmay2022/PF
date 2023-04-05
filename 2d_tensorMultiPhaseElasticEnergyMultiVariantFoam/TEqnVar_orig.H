// scalar sumT0  = 0.0;
// scalar sum2T0 = 0.0;
// scalar sum3T0 = 0.0;
// 
// scalar sumT1  = 0.0;
// scalar sum2T1 = 0.0;
// scalar sum3T1 = 0.0;
// 
// scalar sumT2  = 0.0;
// scalar sum2T2 = 0.0;
// scalar sum3T2 = 0.0;
// 
// scalar sumT3  = 0.0;
// scalar sum2T3 = 0.0;
// scalar sum3T3 = 0.0;



//=================== BEGIN Modified double-well potential =======================//
    volScalarField dwdphi_0  = 18.0*gamma_0*T0*(1.0 -T0)*(1.0 -2.0*T0);
    volScalarField dwdphi_1  = 18.0*gamma_1*T1*(1.0 -T1)*(1.0 -2.0*T1);
    volScalarField dwdphi_2  = 18.0*gamma_2*T2*(1.0 -T2)*(1.0 -2.0*T2);
    volScalarField dwdphi_3  = 18.0*gamma_3*T3*(1.0 -T3)*(1.0 -2.0*T3);
//=================== END Modified double-well potential =========================//


//======================== BEGIN variant evolve equations ===============================================//
    fvScalarMatrix T0Eqn
    (
        sigt*Tau*Epsilon*fvm::ddt(T0) -2.0*gamma_0*Epsilon*sig1*fvm::laplacian(DTs, T0) == -dwdphi_0/(Epsilon)
    );
    
    fvScalarMatrix T1Eqn
    (
        sigt*Tau*Epsilon*fvm::ddt(T1) -2.0*gamma_1*Epsilon*sig1*fvm::laplacian(DTs, T1) == -dwdphi_1/(Epsilon)
    );
    
    fvScalarMatrix T2Eqn
    (
        sigt*Tau*Epsilon*fvm::ddt(T2) -2.0*gamma_2*Epsilon*sig1*fvm::laplacian(DTs, T2) == -dwdphi_2/(Epsilon)
    );
    
    fvScalarMatrix T3Eqn
    (
        sigt*Tau*Epsilon*fvm::ddt(T3) -2.0*gamma_3*Epsilon*sig1*fvm::laplacian(DTs, T3) == -dwdphi_3/(Epsilon)
    );
    
    T0Eqn.solve();
    T1Eqn.solve();
    T2Eqn.solve();
    T3Eqn.solve();
 
//     scalar lambda_phase = 0.0;
    scalar numphases = 4.0;   
    
    lambda_phase = sigt*((fvc::ddt(T0) + fvc::ddt(T1) + fvc::ddt(T2) + fvc::ddt(T3))/numphases);
    Info<< "Min/max lambda_phase:" << min(lambda_phase()).value() << '\t'<< max(lambda_phase()).value() << endl;
//     lambda_phase = ((fvc::ddt(T0) + fvc::ddt(T1) + fvc::ddt(T2) + fvc::ddt(T3))/numphases);
//     lambda_phase /= numphases;

//=============== Substracting lambda_phase from evolution equations ==============================//    
    deltaT0 = sigt*fvc::ddt(T0) - lambda_phase;
    deltaT1 = sigt*fvc::ddt(T1) - lambda_phase;
    deltaT2 = sigt*fvc::ddt(T2) - lambda_phase;
    deltaT3 = sigt*fvc::ddt(T3) - lambda_phase;
//=============== Substracting lambda_phase from evolution equations ==============================//  
    
  //Update the phase-values such that sum upto 1
    T0 +=  deltaT0*runTime.deltaTValue();
    T1 +=  deltaT1*runTime.deltaTValue();
    T2 +=  deltaT2*runTime.deltaTValue();
    T3 +=  deltaT3*runTime.deltaTValue();
  ///////////////////////////////////////////////
   
//==================== BEGIN calculation for volume preservation ==========================//
    scalar lambdaVol_0 = 0.0;
    scalar lambdaVol_1 = 0.0;
    scalar lambdaVol_2 = 0.0;
    
    scalar sumHprime_0 = 0.0;
    scalar sumHprime_1 = 0.0;
    scalar sumHprime_2 = 0.0;
    
    scalar sumRHS0 = 0.0;
    scalar sumRHS1 = 0.0;
    scalar sumRHS2 = 0.0;
	
    forAll(deltaT0, cellI){
      sumRHS0 += deltaT0[cellI];
      sumRHS1 += deltaT1[cellI];
      sumRHS2 += deltaT2[cellI];
      sumHprime_0 += 6.0*T0[cellI]*(1-T0[cellI]);
      sumHprime_1 += 6.0*T1[cellI]*(1-T1[cellI]);
      sumHprime_2 += 6.0*T2[cellI]*(1-T2[cellI]);
    }
    
    lambdaVol_0  = (2.0*sumRHS0 + sumRHS1 + sumRHS2)/sumHprime_0;
//     lambdaVol_0 /= sumHprime_0;
    
    lambdaVol_1  = (sumRHS0 + 2.0*sumRHS1 + sumRHS2)/sumHprime_1;
//     lambdaVol_1 /= sumHprime_1;
    
    lambdaVol_2  = (sumRHS0 + sumRHS1 + 2.0*sumRHS2)/sumHprime_2;
//     lambdaVol_2 /= sumHprime_2;
    Info<< " lambdaVol_0:" << lambdaVol_0 << endl;
    Info<< " lambdaVol_1:" << lambdaVol_1 << endl;
    Info<< " lambdaVol_2:" << lambdaVol_2 << endl;
    
    volScalarField lambdaVol_3 = ((6.0*T0*(1-T0)*lambdaVol_0 + 6.0*T1*(1-T1)*lambdaVol_1 + 6.0*T2*(1-T2)*lambdaVol_2)/numphases);
//     lambdaVol_3 /= numphases;
//============================ END calculation for volume preservation ==================================//

//========== Substracting lambda terms from evolution equation for volume preserve ==========//
    deltaT0 = -(6.0*T0*(1-T0)*lambdaVol_0 -lambdaVol_3);
    deltaT1 = -(6.0*T1*(1-T1)*lambdaVol_1 -lambdaVol_3);
    deltaT2 = -(6.0*T2*(1-T2)*lambdaVol_2 -lambdaVol_3);
    deltaT3 = -lambdaVol_3;
//========== Substracting lambda terms from evolution equation for volume preserve ============//
    
//=========================== phase field update =============================//
//     T0 += (deltaT0 - sigt*fvc::ddt(T0))*runTime.deltaTValue();	//precipitate1
//     T1 += (deltaT1 - sigt*fvc::ddt(T1))*runTime.deltaTValue();	//precipitate2
//     T2 += (deltaT2 - sigt*fvc::ddt(T2))*runTime.deltaTValue();	//precipitate3
//     T3 += (deltaT3 - sigt*fvc::ddt(T3))*runTime.deltaTValue();	//matrix
    T0 +=  deltaT0*runTime.deltaTValue();
    T1 +=  deltaT1*runTime.deltaTValue();
    T2 +=  deltaT2*runTime.deltaTValue();
    T3 +=  deltaT3*runTime.deltaTValue();
//=========================== phase field update =============================//


//================================= END variant evolve equations ====================================//