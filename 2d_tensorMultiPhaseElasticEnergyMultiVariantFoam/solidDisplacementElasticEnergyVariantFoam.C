/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    laplacianFoam

Description
    Solves a simple Laplace equation, e.g. for thermal diffusion in a solid.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fvOptions.H"
#include "simpleControl.H"
#include "Switch.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"

    #include "createTime.H"
    #include "createMesh.H"

    simpleControl simple(mesh);

    #include "createFields.H"

    #include "postProcess.H"
    #include "createControls.H"
   
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nCalculating temperature distribution\n" << endl;

////////////////////Start without elasticity/////////////////
//     while (simple.loop(runTime))
//     {
//         Info<< "Time = " << runTime.timeName() << nl << endl;
//             
//         while (simple.correctNonOrthogonal())
//         {
// // 	  #include "TEqn_2phase.H"
// // 	  #include "TEqnVar.H"
// 	  #include "TEqnVarElast.H"
// 
//         }
// 
// 	runTime.write();
//         Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
//             << "  ClockTime = " << runTime.elapsedClockTime() << " s"
//             << nl << endl;
//     }
//////////////////////End without elasticity/////////////////

/////////////////////Start elasticity///////////////////////////
//     scalar t;
//     for (t=0; t < 10; t++) {
//       #include "TEqnVar_without_elast.H"
//     }
    
    while (runTime.loop())
    {
        Info<< "Iteration: " << runTime.value() << nl << endl;

        #include "readSolidDisplacementFoamControls.H"

        int iCorr = 0;
        scalar initialResidual = 0;
//                 volScalarField T0Old =T0;
// 		volScalarField T1Old =T1;
// 		volScalarField T2Old =T2;
// 		volScalarField T3Old =T3;
//               
            D.correctBoundaryConditions();
            gradD = fvc::grad(D);
          //  gradD = gradD.T(); 
	    #include "TEqnVarElast.H"
            
            do
            {
                 
            {
                
                fvVectorMatrix DEqn
                (
                    fvm::d2dt2(D)
                 ==
                    sig3*fvm::laplacian(2*(mu2*T3*T3*(3-2*T3) + mu1*(1-T3)*(1-T3)*(1+2*T3)) 
		    + lambda2*T3*T3*(3-2*T3) + lambda1*(1-T3)*(1-T3)*(1+2*T3), D, "laplacian(DD,D)")
                  + (sig3/sig2)*divSigmaExp
		  -(sig3)*fvc::div((2*mu1*(1-T3)*(1-T3)*(1+2*T3) + 2*mu2*T3*T3*(3-2*T3))*(T0*T0*(3-2*T0)*cEigenStrain0
                             + T1*T1*(3-2*T1)*cEigenStrain1 +T2*T2*(3-2*T2)*cEigenStrain2)
                             +(lambda2*T3*T3*(3-2*T3) + (1-T3)*(1-T3)*(1+2*T3)*lambda1)*I*tr(T0*T0*(3-2*T0)*cEigenStrain0
                             + T1*T1*(3-2*T1)*cEigenStrain1 +T2*T2*(3-2*T2)*cEigenStrain2))
                );

                /*
                     DEqn -= (sig3)*fvc::div((2*mu1*(1-T3)*(1-T3)*(1+2*T3) + 2*mu2*T3*T3*(3-2*T3))*(T0*T0*(3-2*T0)*cEigenStrain0 
			     + T1*T1*(3-2*T1)*cEigenStrain1 +T2*T2*(3-2*T2)*cEigenStrain2)
			     +(lambda2*T3*T3*(3-2*T3) + (1-T3)*(1-T3)*(1+2*T3)*lambda1)*I*tr(T0*T0*(3-2*T0)*cEigenStrain0 
			     + T1*T1*(3-2*T1)*cEigenStrain1 +T2*T2*(3-2*T2)*cEigenStrain2));
		*/

                initialResidual = DEqn.solve().max().initialResidual();

                if (!compactNormalStress)
                {
                    divSigmaExp = fvc::div(DEqn.flux());
                }
            }

            {
                
                D.correctBoundaryConditions();
                gradD = fvc::grad(D);
		//gradD = gradD.T();
    /*
		scalar eigen_xx=0.0; 
		scalar eigen_yy=0.0;
		scalar eigen_zz=0.0;
		    
                  forAll(strain,cellI)
                  {   
		    eigen_xx  =   T0[cellI]*T0[cellI]*(3-2*T0[cellI])*cEigenStrain0.component(symmTensor::XX).value();
		    eigen_xx +=   T1[cellI]*T1[cellI]*(3-2*T1[cellI])*cEigenStrain1.component(symmTensor::XX).value();
		    eigen_xx +=   T2[cellI]*T2[cellI]*(3-2*T2[cellI])*cEigenStrain2.component(symmTensor::XX).value();
		    
		    eigen_yy  =   T0[cellI]*T0[cellI]*(3-2*T0[cellI])*cEigenStrain0.component(symmTensor::YY).value();
		    eigen_yy +=   T1[cellI]*T1[cellI]*(3-2*T1[cellI])*cEigenStrain1.component(symmTensor::YY).value();
		    eigen_yy +=   T2[cellI]*T2[cellI]*(3-2*T2[cellI])*cEigenStrain2.component(symmTensor::YY).value();
		    
		    eigen_zz  =   T0[cellI]*T0[cellI]*(3-2*T0[cellI])*cEigenStrain0.component(symmTensor::ZZ).value();
		    eigen_zz +=   T1[cellI]*T1[cellI]*(3-2*T1[cellI])*cEigenStrain1.component(symmTensor::ZZ).value();
		    eigen_zz +=   T2[cellI]*T2[cellI]*(3-2*T2[cellI])*cEigenStrain2.component(symmTensor::ZZ).value();
		    
                    strain[cellI].component(symmTensor::XX) = gradD[cellI].component(tensor::XX) - eigen_xx;
                    strain[cellI].component(symmTensor::YY) = gradD[cellI].component(tensor::YY) - eigen_yy;
                    strain[cellI].component(symmTensor::ZZ) = gradD[cellI].component(tensor::ZZ) - eigen_zz;
                    strain[cellI].component(symmTensor::XY) = 0;
                    strain[cellI].component(symmTensor::XZ) = 0;
                    strain[cellI].component(symmTensor::YZ) = 0;
                  }*/
		
   strain=((gradD-((T0*T0*(3-2*T0)*cEigenStrain0)+(T1*T1*(3-2*T1)*cEigenStrain1)+(T2*T2*(3-2*T2)*cEigenStrain2)))
                &&symmTensor(1,0,0,0,0,0))*symmTensor(1,0,0,0,0,0)
          +((gradD-((T0*T0*(3-2*T0)*cEigenStrain0)+(T1*T1*(3-2*T1)*cEigenStrain1)+(T2*T2*(3-2*T2)*cEigenStrain2)))
                &&symmTensor(0,0,0,1,0,0))*symmTensor(0,0,0,1,0,0)
          +((gradD-((T0*T0*(3-2*T0)*cEigenStrain0)+(T1*T1*(3-2*T1)*cEigenStrain1)+(T2*T2*(3-2*T2)*cEigenStrain2)))
                &&symmTensor(0,0,0,0,0,1))*symmTensor(0,0,0,0,0,1);
		
	
                sigmaD = (mu2*T3*T3*(3-2*T3) + mu1*(1-T3)*(1-T3)*(1+2*T3))*twoSymm(gradD) 
		+ (lambda2*T3*T3*(3-2*T3)    + lambda1*(1-T3)*(1-T3)*(1+2*T3))*(I*tr(gradD))
                + (mu2_*T3*T3*(3-2*T3) + mu1_*(1-T3)*(1-T3)*(1+2*T3))*strain;
            
                if (compactNormalStress)
                {
                    divSigmaExp = sig2*fvc::div
                    (
                        sigmaD - (2*mu2*T3*T3*(3-2*T3) + 2*mu1*(1-T3)*(1-T3)*(1+2*T3)
			+lambda2*T3*T3*(3-2*T3) + (1-T3)*(1-T3)*(1+2*T3)*lambda1)*gradD,
                        "div(sigmaD)"
                    );
                }
                else
                {
                    divSigmaExp += sig2*fvc::div(sigmaD);
                }
            }

        } while (initialResidual > convergenceTolerance && ++iCorr < nCorr);
        

        #include "calculateStress.H"

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;

//total energy calculation

    volScalarField consta(0.5*(Sigma && (symm(gradD) -T0*T0*(3-2*T0)*cEigenStrain0 
                          -T1*T1*(3-2*T1)*cEigenStrain1 -T2*T2*(3-2*T2)*cEigenStrain2)));

    volScalarField constc((1-T3)*(0.5*(Sigma && (symm(gradD) -T0*T0*(3-2*T0)*cEigenStrain0
                          -T1*T1*(3-2*T1)*cEigenStrain1 -T2*T2*(3-2*T2)*cEigenStrain2))));

    volVectorField gradT0(fvc::grad(T0));
    volVectorField gradT1(fvc::grad(T1));
    volVectorField gradT2(fvc::grad(T2));
    volVectorField gradT3(fvc::grad(T3));
    
    volScalarField constb(2.0*gamma_0*Epsilon*(magSqr(gradT0)) + 2.0*gamma_1*Epsilon*(magSqr(gradT1)) 
                          + 2.0*gamma_2*Epsilon*(magSqr(gradT2)) + 2.0*gamma_3*Epsilon*(magSqr(gradT3)));
            
    dimensionedScalar elasticEnergyGsum = gSum(consta());
    dimensionedScalar surfaceEnergyGsum = gSum(constb());
    dimensionedScalar totalEnergyGsum   = elasticEnergyGsum + surfaceEnergyGsum;
    dimensionedScalar elasticEnergyPrecipitateGum = gSum(constc());

    Info<< "elasticEnergyGsum: " << (elasticEnergyGsum) << endl;
    Info<< "surfaceEnergyGsum: " << (surfaceEnergyGsum) << endl;
    Info<< "totalEnergyGsum: " << (totalEnergyGsum) << endl;
    Info<< "elasticEnergyPrecipitateGum: " << (elasticEnergyPrecipitateGum) << endl;


// Info<< "Min/max strain:" << min(strain.component(symmTensor::XX)).value() << ' '
// << max(strain.component(symmTensor::XX)).value() << endl;  
// 
// Info<< "Min/max gradD:" << min(gradD.component(tensor::XX)()).value() << ' '
// << max(gradD.component(tensor::XX)()).value() << endl;

//for total energy calculation//
//     scalar sumTotalEnergy = 0.0;
/*    dimensionedScalar totalEnergy = 0.0;
//     scalar consta = 0.0;
//     scalar constb = 0.0;
    volScalarField consta(Sigma && (symm(gradD)-T*T*(3-2*T)*cEigenStrain));//symm(fvc::grad(D));
//     volScalarField constb(Sigma && cEigenStrain);
    volVectorField gradT(fvc::grad(T));
    forAll(consta, cellI) {
	totalEnergy += 0.5*consta[cellI] + 2.0*Gamma*Epsilon*pow(mag(gradT[cellI]),2);
     }*/   
	    
// Info<< "totalEnergy:" << totalEnergy << endl;
	    
// Info<< "Min/max T:" << min(T()).value() << ' '
//     << max(T()).value() << endl;
    }

///////////////////////End elasticity//////////////////////////////////////////




    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
