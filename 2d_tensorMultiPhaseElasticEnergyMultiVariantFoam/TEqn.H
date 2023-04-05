/////////////////////single precipitate///////////////////////////////////
scalar sumT = 0.0;
scalar sum2T = 0.0;
scalar sum3T = 0.0;
   
forAll(T,cellI) {
    strain[cellI].component(symmTensor::XX) = gradD[cellI].component(tensor::XX) - T[cellI]*T[cellI]*(3-2*T[cellI])*cEigenStrain.component(symmTensor::XX).value();
    strain[cellI].component(symmTensor::YY) = gradD[cellI].component(tensor::YY) - T[cellI]*T[cellI]*(3-2*T[cellI])*cEigenStrain.component(symmTensor::YY).value();
    strain[cellI].component(symmTensor::ZZ) = gradD[cellI].component(tensor::ZZ) - T[cellI]*T[cellI]*(3-2*T[cellI])*cEigenStrain.component(symmTensor::ZZ).value();
    strain[cellI].component(symmTensor::XY) = 0;
    strain[cellI].component(symmTensor::XZ) = 0;
    strain[cellI].component(symmTensor::YZ) = 0;
}
   
   
   Sigma = (2*(mu1*T*T*(3-2*T) + mu2*(1-T)*(1-T)*(1+2*T))*(symm(fvc::grad(D)) - T*T*(3-2*T)*cEigenStrain) 
    + (lambda1*T*T*(3-2*T) + lambda2*(1-T)*(1-T)*(1+2*T))*(I*tr(fvc::grad(D) - T*T*(3-2*T)*cEigenStrain)))
    + (mu1_*T*T*(3-2*T) + mu2_*(1-T)*(1-T)*(1+2*T))*strain;


   deltaSigmaD = ((mu1-mu2)*twoSymm(fvc::grad(D)) + (lambda1-lambda2)*(I*tr(fvc::grad(D))) 
    - 2*(mu1-mu2)*T*T*(3-2*T)*cEigenStrain - (lambda1-lambda2)*(I*tr(T*T*(3-2*T)*cEigenStrain)))
    + (mu1_-mu2_)*strain;
    
 
    fvScalarMatrix TEqn
    (
        Tau*Epsilon*fvm::ddt(T)- 2.0*Gamma*Epsilon*fvm::laplacian(DTs, T)== -(18.0*Gamma/(sig*Epsilon))*(T)*(1-T)*(1-2.0*T)
    );
    
    TEqn -= (1/sig)*6.0*T*(1-T)*(0.5*(deltaSigmaD && (symm(fvc::grad(D))-T*T*(3-2*T)*cEigenStrain))-(Sigma && cEigenStrain));

    TEqn.solve();

forAll(T, cellI) {
    sumT += T[cellI];
    sum2T+= pow(T[cellI],2);
    sum3T += T[cellI]-TOld[cellI];
}

scalar lambdas=sum3T/(6*(sumT-sum2T));

forAll(T, cellI) {
  T[cellI] += - 6*lambdas*T[cellI]*(1-T[cellI]);
}
/////////////////////single precipitate///////////////////////////////////


