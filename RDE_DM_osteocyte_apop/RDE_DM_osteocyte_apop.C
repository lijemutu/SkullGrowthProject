/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Application
    RDE_DyM

Description
    Reaction-Diffusion-Strain solver with dynamic mesh technique for cranial vault growth.

Author
    Chanyoung Lee

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "constitutiveModel.H"
#include "solidInterface.H"
#include "volPointInterpolation.H"
#include "pointPatchInterpolation.H"
#include "primitivePatchInterpolation.H"
#include "pointFields.H"
#include "twoDPointCorrector.H"
#include "leastSquaresVolPointInterpolation.H"
#include "processorFvPatchFields.H"
#include "transformGeometricField.H"
#include "symmetryPolyPatch.H"
#include "dynamicFvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
# include "setRootCase.H"
# include "createTime.H"
//# include "createMesh.H"
#   include "createDynamicFvMesh.H"
# include "createFields.H"
# include "readDivDSigmaExpMethod.H"
# include "readDivDSigmaLargeStrainExpMethod.H"
# include "readMoveMeshMethod.H"
# include "createSolidInterfaceNonLin.H"
# include "findGlobalFaceZones.H"

// Read check frequency
    label checkFrequency = 1;
    args.optionReadIfPresent("checkFrequency", checkFrequency);
    
//* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

  Info << "\nStarting time loop\n" << endl;
  

  for (runTime++; !runTime.end(); runTime++)
    {
      Info<< "Time = " << runTime.timeName() << nl << endl;
      


#     include "readSolidMechanicsControls.H"

rho = (rho_m-rho_b)*(pow(oS,m)/(pow(o,m)+pow(oS,m))) + rho_b;
rho_v = rho;
E_v = (Em-Eb)*(pow(oS,m)/(pow(o,m)+pow(oS,m))) + Eb;
nu_v = (nu_m-nu_b)*(pow(oS,m)/(pow(o,m)+pow(oS,m))) + nu_b;

mu = E_v/(2*(1+nu_v));
lambda = nu_v*E_v/((1+nu_v)*(1-2*nu_v));
muf= fvc::interpolate(mu);
lambdaf = fvc::interpolate(lambda); 





      int iCorr = 0;
      lduSolverPerformance solverPerf;
      scalar initialResidual = 1.0;
      scalar relativeResidual = 1.0;
      lduMatrix::debug = 0;

      
      
      do
      {
          DU.storePrevIter();

#         include "calculateDivDSigmaExp.H"
#         include "calculateDivDSigmaLargeStrainExp.H"

          //- Updated lagrangian momentum equation
          fvVectorMatrix DUEqn
              (
                  fvm::d2dt2(rho,DU)
                  ==
                  fvm::laplacian(2*muf + lambdaf, DU, "laplacian(DDU,DU)")
                  + divDSigmaExp
                  + divDSigmaLargeStrainExp
                  );



          if (solidInterfaceCorr)
          {
              solidInterfacePtr->correct(DUEqn);
          }

          solverPerf = DUEqn.solve();

          if (iCorr == 0)
          {
              initialResidual = solverPerf.initialResidual();
          }
	  
          DU.relax();
	  

          gradDU = fvc::grad(DU);
	  


#         include "calculateDEpsilonDSigma.H"
#         include "calculateRelativeResidual.H"



          Info << "\tTime " << runTime.value()
               << ", Corrector " << iCorr
               << ", Solving for " << DU.name()
               << " using " << solverPerf.solverName()
               << ", res = " << solverPerf.initialResidual()
               << ", rel res = " << relativeResidual
               << ", inner iters " << solverPerf.nIterations() << endl;
      }
      while
        (
            //solverPerf.initialResidual() > convergenceTolerance
            relativeResidual > convergenceTolerance
            && ++iCorr < nCorr
            );

      Info << nl << "Time " << runTime.value() << ", Solving for " << DU.name()
           << ", Initial residual = " << initialResidual
           << ", Final residual = " << solverPerf.initialResidual()
           << ", No outer iterations " << iCorr << endl;


#     include "moveMesh.H"

#     include "rotateFields.H"

#     include "writeFields.H"


      Info<< nl << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
          << "  ClockTime = " << runTime.elapsedClockTime() << " s"
          << endl;



volTensorField F = I + gradDU;

volTensorField Finv = inv(F);

DEV = det(F)-1;
REV = DEV/(runTime.time().deltaT());
EV += DEV;

volVectorField v = DU/(runTime.time().deltaT()); 
surfaceScalarField phi = linearInterpolate(v) & mesh.Sf();

#       include "readSIMPLEControls.H"

	Da_v = (pow(oS,m)/(pow(o,m)+pow(oS,m)))*Da;
	Dh_v = (pow(oS,m)/(pow(o,m)+pow(oS,m)))*Dh;

	      
        for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
        {
             solve
            (
                fvm::ddt(a) 
		-a*fvc::div(phi)

		-(pow(oS,m)/(pow(o,m)+pow(oS,m)))*((RE0+REV))*alpha_a
	
		-(pow(oS,m)/(pow(o,m)+pow(oS,m)))*((RE0+REV))*alpha_o*o
		+(pow(oS,m)/(pow(o,m)+pow(oS,m)))*beta_a*a
		
		-(pow(oS,m)/(pow(o,m)+pow(oS,m)))*gamma_a*a*a/h
		
		- fvm::laplacian(Da_v, a)
				
            );
	    
	    solve
            (
                fvm::ddt(h) 
		
		-h*fvc::div(phi)
		-(pow(oS,m)/(pow(o,m)+pow(oS,m)))*alpha_h
		+(pow(oS,m)/(pow(o,m)+pow(oS,m)))*beta_h*h
		
		
		- (pow(oS,m)/(pow(o,m)+pow(oS,m)))*gamma_h*a*a
		- fvm::laplacian(Dh_v, h)
		
            );
	
            solve
            (
                fvm::ddt(o)
		-o*fvc::div(phi)
               
		-eta*(pow(oS,m)/(pow(o,m)+pow(oS,m)))*(pow((a*a/h),m)/(pow((a*a/h),m)+pow(aT,m)))*(pow(EV,m)/(pow(EV,m)+pow(ET,m)))*(pow(b_s,m)/(pow(b_s,m)+pow(b,m)))+etaap*a*(pow(aT,m)/(pow(aT,m)+pow((a*a/h),m)))+etaapb*a*(pow(86400,m)/(pow(runTime.time().value()-tau_1,m)+pow(86400,m)))

		
            ); 


            solve
            (
                fvm::ddt(ocy)
		-ocy*fvc::div(phi)
               
		-eta_ocy*(pow(o,m)/(pow(oS,m)+pow(o,m)))*(1+D_T*EV)
		+eta_apocy*(pow(aT,m)/(pow(aT,m)+pow((a*a/h),m)))*(pow(runTime.time().value(),m)/(pow(tau_2,m)+pow(runTime.time().value(),m)))
		
            ); 


            solve
            (
                fvm::ddt(b)
		-b*fvc::div(phi)
               
		-eta_b*(pow(o,m)/(pow(oS,m)+pow(o,m)))
            );



	   
        }
	

oGrad = mag(fvc::grad(o));
    
	     

        runTime.write();
	mesh.update();


    }

  Info<< "End\n" << endl;

  return(0);
}

// ************************************************************************* //
