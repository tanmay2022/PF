#include "readMechanicalProperties.H"
#include "readThermalProperties.H"

Info<< "Reading field D\n" << endl;
IOdictionary constProperties
    (
        IOobject
        (
            "constProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );
    
    dimensionedScalar gamma_0
    (
        constProperties.lookup("gamma_0")
    );    
    dimensionedScalar gamma_1
    (
        constProperties.lookup("gamma_1")
    );    
    dimensionedScalar gamma_2
    (
        constProperties.lookup("gamma_2")
    );    
    dimensionedScalar gamma_3
    (
        constProperties.lookup("gamma_3")
    );        
    dimensionedScalar sigt
    (
	constProperties.lookup("sigt")
    );
    dimensionedScalar Tau
    (
        constProperties.lookup("Tau")
    );
    dimensionedScalar Epsilon
    (
        constProperties.lookup("Epsilon")
    );
    
//=====================for three variants============================//
    dimensionedScalar DTs
    (
        constProperties.lookup("DTs")
    );
   
    dimensionedScalar sig1
    (
        constProperties.lookup("sig1")
    );
  
    dimensionedScalar lambdas			//for phase balance
    (
        constProperties.lookup("lambdas")
    );
    
    dimensionedSymmTensor cEigenStrain0
    (
        constProperties.lookup("cEigenStrain0")
    );  
    dimensionedSymmTensor cEigenStrain1
    (
        constProperties.lookup("cEigenStrain1")
    );  
    dimensionedSymmTensor cEigenStrain2
    (
        constProperties.lookup("cEigenStrain2")
    );      
/*    dimensionedSymmTensor cEigenStrain3
    (
        constProperties.lookup("cEigenStrain3")
    ); */     

    dimensionedScalar mu1
    (
        constProperties.lookup("mu1")
    );
    dimensionedScalar mu2
    (
        constProperties.lookup("mu2")
    );
    dimensionedScalar mu1_
    (
        constProperties.lookup("mu1_")
    );    
    dimensionedScalar mu2_
    (
        constProperties.lookup("mu2_")
    );    
    dimensionedScalar lambda1
    (
        constProperties.lookup("lambda1")
    );
    dimensionedScalar lambda2
    (
        constProperties.lookup("lambda2")
    );
//     dimensionedScalar sig1
//     (
//         constProperties.lookup("sig1")
//     );
    dimensionedScalar sig2
    (
        constProperties.lookup("sig2")
    );
    dimensionedScalar sig3
    (
        constProperties.lookup("sig3")
    );    
    
volVectorField D
(
    IOobject
    (
        "D",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volTensorField gradD
(
    IOobject
    (
        "gradD",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
//     mesh
    fvc::grad(D)
);

volSymmTensorField strain
(
    IOobject
    (
        "strain",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
    //volSymmTensorField("strain", dimensionSet(0,0,0,0,0,0,0), Foam::symmTensor(0,0,0,0,0,0))
   // Foam::symmTensor(0,0,0,0,0,0)
);

// volScalarField totalEnergy
// (
//     IOobject
//     (
//         "totalEnergy",
//         runTime.timeName(),
//         mesh,
//         IOobject::MUST_READ,
//         IOobject::AUTO_WRITE
//     ),
//     mesh
// );
volSymmTensorField sigma
            (

                IOobject
                (
                    "sigma",
                    runTime.timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh
            );

//===================for single precipitate========================//
// volScalarField T
//         (
//             IOobject
//             (
//                 "T",
//                 runTime.timeName(),
//                 mesh,
//                 IOobject::MUST_READ,
//                 IOobject::AUTO_WRITE
//             ),
//             mesh
//         );
//=====================for three variants============================//
volScalarField T0
        (
            IOobject
            (
                "T0",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh
        );
	
volScalarField T1
        (
            IOobject
            (
                "T1",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh
        );	
	
volScalarField T2
        (
            IOobject
            (
                "T2",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh
        );	

volScalarField T3
        (
            IOobject
            (
                "T3",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh
        );	

// volScalarField T0Old
//     (
//         IOobject
//             (
//                 "T0Old",
//                 runTime.timeName(),
//                 mesh,
//                 IOobject::NO_READ,
//                 IOobject::NO_WRITE
//             ),
//             T0
//     );
// 
// volScalarField T1Old
//     (
//         IOobject
//             (
//                 "T1Old",
//                 runTime.timeName(),
//                 mesh,
//                 IOobject::NO_READ,
//                 IOobject::NO_WRITE
//             ),
//             T1
//     ); 
//     
// volScalarField T2Old
//     (
//         IOobject
//             (
//                 "T2Old",
//                 runTime.timeName(),
//                 mesh,
//                 IOobject::NO_READ,
//                 IOobject::NO_WRITE
//             ),
//             T2
//     );     
//     
// volScalarField T3Old
//     (
//         IOobject
//             (
//                 "T3Old",
//                 runTime.timeName(),
//                 mesh,
//                 IOobject::NO_READ,
//                 IOobject::NO_WRITE
//             ),
//             T3
//     );    
	
volScalarField lambda_phase
        (
            IOobject
            (
                "lambda_phase",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            sigt*((fvc::ddt(T0) + fvc::ddt(T1) + fvc::ddt(T2) + fvc::ddt(T3))/4.0)
        );
	
/*volScalarField lambdaVol_3
        (
            IOobject
            (
                "lambdaVol_3",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            ((6.0*T0*(1-T0)*lambdaVol_0 + 6.0*T1*(1-T1)*lambdaVol_1 + 6.0*T2*(1-T2)*lambdaVol_2)/4.0)
        );*/	
	
volScalarField deltaT0
        (
            IOobject
            (
                "deltaT0",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            (sigt*fvc::ddt(T0) - lambda_phase)
        );
	
volScalarField deltaT1
        (
            IOobject
            (
                "deltaT1",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            (sigt*fvc::ddt(T1) - lambda_phase)
        );	
	
volScalarField deltaT2
        (
            IOobject
            (
                "deltaT2",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            (sigt*fvc::ddt(T2) - lambda_phase)
        );	

volScalarField deltaT3
        (
            IOobject
            (
                "deltaT3",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            (sigt*fvc::ddt(T3) - lambda_phase)
        ); 

// volScalarField lambdaVol_0
//         (
//             IOobject
//             (
//                 "lambdaVol_0",
//                 runTime.timeName(),
//                 mesh,
//                 IOobject::MUST_READ,
//                 IOobject::NO_WRITE
//             ),
//             mesh
//         );
// 	
// volScalarField lambdaVol_1
//         (
//             IOobject
//             (
//                 "lambdaVol_1",
//                 runTime.timeName(),
//                 mesh,
//                 IOobject::MUST_READ,
//                 IOobject::NO_WRITE
//             ),
//             mesh
//         );
// 	
// volScalarField lambdaVol_2
//         (
//             IOobject
//             (
//                 "lambdaVol_2",
//                 runTime.timeName(),
//                 mesh,
//                 IOobject::MUST_READ,
//                 IOobject::NO_WRITE
//             ),
//             mesh
//         );
// 	
// volScalarField lambdaVol_3
//         (
//             IOobject
//             (
//                 "lambdaVol_3",
//                 runTime.timeName(),
//                 mesh,
//                 IOobject::MUST_READ,
//                 IOobject::NO_WRITE
//             ),
//             mesh
//         );	
// 	

//=====================for three variants============================//
	

// volScalarField TOld
//     (
//         IOobject
//             (
//                 "TOld",
//                 runTime.timeName(),
//                 mesh,
//                 IOobject::MUST_READ,
//                 IOobject::AUTO_WRITE
//             ),
//             mesh
//     );




    Info<< "Reading diffusivity DT\n" << endl;

    
//=====================for three variants============================//
//     dimensionedScalar lambdaVol_0		//for volume constraint 
//     (
//         constProperties.lookup("lambdaVol_0")
//     );
// 
//     dimensionedScalar lambdaVol_1		//for volume constraint
//     (
//         constProperties.lookup("lambdaVol_1")
//     );    
//     
//     dimensionedScalar lambdaVol_2		//for volume constraint
//     (
//         constProperties.lookup("lambdaVol_2")
//     );    
//     
//     dimensionedScalar lambdaVol_3		//for volume constraint
//     (
//         constProperties.lookup("lambdaVol_3")
//     );    
//     dimensionedScalar lambda_phase
//     (
//         constProperties.lookup("lambda_phase")
//     );  
    
//     dimensionedSymmTensor cEigenStrain0
//     (
//         constProperties.lookup("cEigenStrain0")
//     );      
// 
//     dimensionedSymmTensor cEigenStrain1
//     (
//         constProperties.lookup("cEigenStrain1")
//     );
// 
//     dimensionedSymmTensor cEigenStrain2
//     (
//         constProperties.lookup("cEigenStrain2")
//     );
//     dimensionedSymmTensor cEigenStrain3
//     (
//         constProperties.lookup("cEigenStrain3")
//     );    
    
//     dimensionedSymmTensor cEigenStrain
//     (
//         constProperties.lookup("cEigenStrain")
//     );     

  

volSymmTensorField sigmaD
(
    IOobject
    (
        "sigmaD",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
//     (mu1*T+mu2*(1-T))*twoSymm(fvc::grad(D)) + (lambda1*T+lambda2*(1-T))*(I*tr(fvc::grad(D)))
    (mu2*T3*T3*(3-2*T3) + mu1*(1-T3)*(1-T3)*(1+2*T3))*twoSymm(gradD) + (lambda2*T3*T3*(3-2*T3)    
    + lambda1*(1-T3)*(1-T3)*(1+2*T3))*(I*tr(gradD))
);

volSymmTensorField deltaSigmaD
(
    IOobject
    (
        "deltaSigmaD",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
//      ((mu1-mu2)*twoSymm(fvc::grad(D)) + (lambda1-lambda2)*(I*tr(fvc::grad(D))) - 2*(mu1-mu2)*T*cEigenStrain - (lambda1-lambda2)*(I*tr(T*cEigenStrain)))
	((mu1)*twoSymm(fvc::grad(D)) + (lambda1)*(I*tr(fvc::grad(D))) 
	- 2*(mu1)*(T0*T0*(3-2*T0)*cEigenStrain0 + T1*T1*(3-2*T1)*cEigenStrain1 +T2*T2*(3-2*T2)*cEigenStrain2) 
	- (lambda1)*(I*tr((T0*T0*(3-2*T0)*cEigenStrain0 + T1*T1*(3-2*T1)*cEigenStrain1 +T2*T2*(3-2*T2)*cEigenStrain2))))
	+ (mu1_)*strain
);
// 
volSymmTensorField Sigma
(
    IOobject
    (
        "Sigma",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
//     (2*(mu1*T+mu2*(1-T))*(0.5*twoSymm(fvc::grad(D)) - T*cEigenStrain) + (lambda1*T+lambda2*(1-T))*(I*tr(fvc::grad(D) - T*cEigenStrain)))
    (2*(mu2*T3*T3*(3-2*T3) + mu1*(1-T3)*(1-T3)*(1+2*T3))*(symm(fvc::grad(D)) 
	- (T0*T0*(3-2*T0)*cEigenStrain0 + T1*T1*(3-2*T1)*cEigenStrain1 +T2*T2*(3-2*T2)*cEigenStrain2)) 
    + (lambda2*T3*T3*(3-2*T3) + lambda1*(1-T3)*(1-T3)*(1+2*T3))*(I*tr(fvc::grad(D) 
	- (T0*T0*(3-2*T0)*cEigenStrain0 + T1*T1*(3-2*T1)*cEigenStrain1 +T2*T2*(3-2*T2)*cEigenStrain2)))) 
    + (mu2_*T3*T3*(3-2*T3) + mu1_*(1-T3)*(1-T3)*(1+2*T3))*strain
);


Info<< "Calculating explicit part of div(sigma) divSigmaExp\n" << endl;
volVectorField divSigmaExp
(
    IOobject
    (
        "divSigmaExp",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    sig2*fvc::div(sigmaD)
);

// if (compactNormalStress)
// {
// //     divSigmaExp -= sig2*fvc::laplacian(2*(mu1*T+mu2*(1-T)) + lambda1*T+lambda2*(1-T), D, "laplacian(DD,D)");
//     
// }
// else
// {
// //     divSigmaExp -= sig2*fvc::div((2*(mu1*T+mu2*(1-T)) + lambda1*T+lambda2*(1-T))*fvc::grad(D), "div(sigmaD)");
// }

mesh.setFluxRequired(D.name());
