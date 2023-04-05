    if (runTime.writeTime())
    {
/*            
            volSymmTensorField sigma
            (
            
                IOobject
                (
                    "sigma",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                ((2*(mu1*T*T*(3-2*T) + mu2*(1-T)*(1-T)*(1+2*T))*(symm(fvc::grad(D)) - T*T*(3-2*T)*cEigenStrain)
              + (lambda1*T*T*(3-2*T) + lambda2*(1-T)*(1-T)*(1+2*T))*(I*tr(fvc::grad(D) - T*T*(3-2*T)*cEigenStrain)))
              + (mu1_*T*T*(3-2*T) + mu2_*(1-T)*(1-T)*(1+2*T))*strain)
            );
            
//             sigma = sigma - (lambda1*T+lambda2*(1-T))*rho*I*tr(T*cEigenStrain)-(rho*(2*(mu1*T+mu2*(1-T))*T*cEigenStrain));
*/            
           /* sigma = ((2*(mu1*T*T*(3-2*T) + mu2*(1-T)*(1-T)*(1+2*T))*(symm(fvc::grad(D)) - T*T*(3-2*T)*cEigenStrain) 
            + (lambda1*T*T*(3-2*T) + lambda2*(1-T)*(1-T)*(1+2*T))*(I*tr(fvc::grad(D) - T*T*(3-2*T)*cEigenStrain)))
            + (mu1_*T*T*(3-2*T) + mu2_*(1-T)*(1-T)*(1+2*T))*strain);*/
/*     
        volScalarField sigmaEq
        (
            IOobject
            (
                "sigmaEq",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            sqrt((3.0/2.0)*magSqr(dev(sigma)))
        );

        Info<< "Max sigma = " << max(sigma).value()
            << endl;
        Info<< "Max sigmaEq = " << max(sigmaEq).value()
            << endl;
*/
        runTime.write();

    }
