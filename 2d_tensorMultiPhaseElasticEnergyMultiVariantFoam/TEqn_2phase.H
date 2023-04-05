/////////////////////single precipitate///////////////////////////////////
scalar sumT = 0.0;
scalar sum2T = 0.0;
scalar sum3T = 0.0;
scalar sumT3 = 0.0;
scalar sum2T3 = 0.0;
scalar sum3T3 = 0.0;
    
 
    fvScalarMatrix T0Eqn
    (
        sigt*Tau*Epsilon*fvm::ddt(T0) -2.0*gamma_0*Epsilon*sig1*fvm::laplacian(DTs, T0)== -(18.0*gamma_0*T0*(1.0 -T0)*(1.0 -2.0*T0))/Epsilon
    );
    T0Eqn.solve();
    fvScalarMatrix T3Eqn
    (
        sigt*Tau*Epsilon*fvm::ddt(T3) -2.0*gamma_3*Epsilon*sig1*fvm::laplacian(DTs, T3)== -(18.0*gamma_3*T3*(1.0 -T3)*(1.0 -2.0*T3))/Epsilon
    );   
    T3Eqn.solve();

    scalar numphases = 2.0;   
    
    lambda_phase = ((T0-T0Old)+(T3-T3Old))/numphases;
    
    T0 += T0 -T0Old - lambda_phase;
    T3 += T3 -T3Old - lambda_phase;
    
    
    forAll(T0, cellI) {
	sumT += T0[cellI];
	sum2T+= pow(T0[cellI],2);
	sum3T += T0[cellI]-T0Old[cellI];
    }
    
//      forAll(T3, cellI) {
// 	sumT3 += T3[cellI];
// 	sum2T3+= pow(T3[cellI],2);
// 	sum3T3 += T3[cellI]-T3Old[cellI];
//     }

    scalar lambdaVol_0=sum3T/(6*(sumT-sum2T));
   // scalar lambdaVol_3=sum3T3/(6*(sumT3-sum2T3));

//     scalar lambdaVol_3 = 0;
//     scalar lambdaVol_3 = ((6.0*T0*(1-T0)*lambdaVol_0)/numphases);
    
//     forAll(T3, cellI){
//       lambdaVol_3 = (6.0*T0[cellI]*(1-T0[cellI])*lambdaVol_0)/numphases;
//     }
    
    Info<< "lambdaVol_0:" << lambdaVol_0 << endl;
//     Info<< "lambdaVol_3:" << lambdaVol_3 << endl;    
    
    forAll(T0, cellI) {
      T0[cellI] += - 6*lambdaVol_0*T0[cellI]*(1-T0[cellI]);
      T3[cellI] += (6.0*T0[cellI]*(1-T0[cellI])*lambdaVol_0);
    }
    
/////////////////////single precipitate///////////////////////////////////


