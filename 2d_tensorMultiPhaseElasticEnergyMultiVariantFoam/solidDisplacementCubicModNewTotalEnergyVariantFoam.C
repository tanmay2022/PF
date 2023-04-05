/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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
    solidDisplacementFoam

Description
    Transient segregated finite-volume solver of linear-elastic,
    small-strain deformation of a solid body, with optional thermal
    diffusion and thermal stresses.

    Simple linear elasticity structural analysis code.
    Solves for the displacement vector field D, also generating the
    stress tensor field sigma.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "postProcess.H"

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControls.H"
    #include "createFields.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nCalculating displacement field\n" << endl;
  

    while (runTime.loop())
    {
        Info<< "Iteration: " << runTime.value() << nl << endl;
	#include "TEqnVar.H"
	/*
        #include "readSolidDisplacementFoamControls.H"

        int iCorr = 0;
        scalar initialResidual = 0;
        
        
        
        forAll(T, cellI) 
            {
                TOld[cellI]=T[cellI];
            }   
            
            gradD = fvc::grad(D);
            
            
            do
            {
                 
            {
                
                fvVectorMatrix DEqn
                (
                    fvm::d2dt2(D)
                 ==
                    sig1*fvm::laplacian(2*(mu1*T*T*(3-2*T) + mu2*(1-T)*(1-T)*(1+2*T)) 
		    + lambda1*T*T*(3-2*T) + lambda2*(1-T)*(1-T)*(1+2*T), D, "laplacian(DD,D)")
                  + (sig1/sig2)*divSigmaExp
                );

                
                     DEqn -= (sig1)*fvc::div((2*mu1*T*T*(3-2*T) + 2*mu2*(1-T)*(1-T)*(1+2*T))*T*T*(3-2*T)*cEigenStrain
                     +(lambda1*T*T*(3-2*T) + (1-T)*(1-T)*(1+2*T)*lambda2)*I*tr(T*T*(3-2*T)*cEigenStrain));
                 
		     
                initialResidual = DEqn.solve().max().initialResidual();

                if (!compactNormalStress)
                {
                    divSigmaExp = fvc::div(DEqn.flux());
                }
            }

            {
                
                gradD = fvc::grad(D);

                  forAll(strain,cellI)
                  {
                    strain[cellI].component(symmTensor::XX) = gradD[cellI].component(tensor::XX) - T[cellI]*T[cellI]*(3-2*T[cellI])*cEigenStrain.component(symmTensor::XX).value();
                    strain[cellI].component(symmTensor::YY) = gradD[cellI].component(tensor::YY) - T[cellI]*T[cellI]*(3-2*T[cellI])*cEigenStrain.component(symmTensor::YY).value();
                    strain[cellI].component(symmTensor::ZZ) = gradD[cellI].component(tensor::ZZ) - T[cellI]*T[cellI]*(3-2*T[cellI])*cEigenStrain.component(symmTensor::ZZ).value();
                    strain[cellI].component(symmTensor::XY) = 0;
                    strain[cellI].component(symmTensor::XZ) = 0;
                    strain[cellI].component(symmTensor::YZ) = 0;
                  }

                //  DEqn += (sig1)*fvc::div((mu1_*T*T*(3-2*T) + mu2_*(1-T)*(1-T)*(1+2*T))*strain);
                //

                sigmaD = (mu1*T*T*(3-2*T) + mu2*(1-T)*(1-T)*(1+2*T))*twoSymm(gradD) 
		+ (lambda1*T*T*(3-2*T)    + lambda2*(1-T)*(1-T)*(1+2*T))*(I*tr(gradD))
                + (mu1_*T*T*(3-2*T) + mu2_*(1-T)*(1-T)*(1+2*T))*strain;
            
                if (compactNormalStress)
                {
                    divSigmaExp = sig2*fvc::div
                    (
                        sigmaD - (2*mu1*T*T*(3-2*T) + 2*mu2*(1-T)*(1-T)*(1+2*T)
			+lambda1*T*T*(3-2*T) + (1-T)*(1-T)*(1+2*T)*lambda2)*gradD,
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


    dimensionedScalar totalEnergy = 0.0;

    volScalarField consta(Sigma && (symm(gradD)-T*T*(3-2*T)*cEigenStrain));//symm(fvc::grad(D));

    volVectorField gradT(fvc::grad(T));
    forAll(consta, cellI) {
	totalEnergy += 0.5*consta[cellI] + 2.0*Gamma*Epsilon*pow(mag(gradT[cellI]),2);
     }   
	    
Info<< "totalEnergy:" << totalEnergy << endl;


	    
Info<< "Min/max T:" << min(T()).value() << ' '
    << max(T()).value() << endl;
    }
*/
       runTime.write();
    }
   
    Info<< "End\n" << endl;   
    
    return 0;
}


// ************************************************************************* //
