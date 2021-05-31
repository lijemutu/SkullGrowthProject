 

    Info<< "Reading incremental displacement field DU\n" << endl;
    volVectorField DU
    (
        IOobject
        (
            "DU",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    volTensorField gradDU = fvc::grad(DU);

    Info<< "Reading accumulated displacement field U\n" << endl;
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector("zero", dimLength, vector::zero)
    );
    

    Info << "Reading accumulated strain field epsilon\n" << endl;
    volSymmTensorField epsilon
    (
        IOobject
        (
            "epsilon",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("zero", dimless, symmTensor::zero)
    );



    volSymmTensorField DEpsilon
    (
        IOobject
        (
            "DEpsilon",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("zero", dimless, symmTensor::zero)
    );

    Info << "Reading accumulated stress field sigma\n" << endl;
    volSymmTensorField sigma
    (
        IOobject
        (
            "sigma",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("zero", dimForce/dimArea, symmTensor::zero)
    );


    Info << "Reading incremental stress field DSigma\n" << endl;
    volSymmTensorField DSigma
    (
        IOobject
        (
            "DSigma",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("zero", dimForce/dimArea, symmTensor::zero)
    );

    //- explicit terms in the momentum equation
    volVectorField divDSigmaExp
	(
	 IOobject
	 (
	  "divDSigmaExp",
	  runTime.timeName(),
	  mesh,
	  IOobject::NO_READ,
	  IOobject::NO_WRITE
	  ),
	 mesh,
 	 dimensionedVector("zero", dimensionSet(1, -2, -2, 0, 0, 0, 0), vector::zero)
	 );

     volVectorField divDSigmaLargeStrainExp
	(
	 IOobject
	 (
	  "divDSigmaLargeStrainExp",
	  runTime.timeName(),
	  mesh,
	  IOobject::NO_READ,
	  IOobject::NO_WRITE
	  ),
	 mesh,
 	 dimensionedVector("zero", dimensionSet(1, -2, -2, 0, 0, 0, 0), vector::zero)
	 );

    constitutiveModel rheology(sigma, DU);

    volScalarField rho = rheology.rho();

 
volScalarField rho_v
    (
        IOobject
        (
            "rho_v",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );


	
    volScalarField mu = rheology.mu();
    volScalarField lambda = rheology.lambda();
    surfaceScalarField muf = fvc::interpolate(rheology.mu());
    surfaceScalarField lambdaf = fvc::interpolate(rheology.lambda());

    surfaceVectorField n = mesh.Sf()/mesh.magSf();





    volScalarField E_v
    (
        IOobject
        (
            "E_v",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );


    volScalarField nu_v
    (
        IOobject
        (
            "nu_v",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
    
    

    Info<< "Reading field a\n" << endl;

    volScalarField a
    (
        IOobject
        (
            "a",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
 



  Info<< "Reading field h\n" << endl;

    volScalarField h
    (
        IOobject
        (
            "h",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );


Info<< "Reading field o\n" << endl;

    volScalarField o
    (
        IOobject
        (
            "o",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
    
volScalarField oGrad
(IOobject
("oGrad",
runTime.timeName(),
mesh,
IOobject::READ_IF_PRESENT,
IOobject::AUTO_WRITE
),
mag(fvc::grad(o))
);


 volScalarField Da_v
    (
        IOobject
        (
            "Da_v",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    
    );
    

volScalarField Dh_v
    (
        IOobject
        (
            "Dh_v",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    
    );


volScalarField b
    (
        IOobject
        (
            "b",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    
    );

volScalarField ocy
    (
        IOobject
        (
            "ocy",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    
    );
    
    



    Info<< "Reading transportProperties\n" << endl;

    IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

Info<< "Reading alpha_a\n" << endl;

    dimensionedScalar alpha_a
    (
        transportProperties.lookup("alpha_a")
    );

Info<< "Reading alpha_h\n" << endl;

    dimensionedScalar alpha_h
    (
        transportProperties.lookup("alpha_h")
    );

Info<< "Reading alpha_o\n" << endl;

    dimensionedScalar alpha_o
    (
        transportProperties.lookup("alpha_o")
    );


Info<< "Reading beta_a\n" << endl;

    dimensionedScalar beta_a
    (
        transportProperties.lookup("beta_a")
    );

Info<< "Reading beta_h\n" << endl;

    dimensionedScalar beta_h
    (
        transportProperties.lookup("beta_h")
    );




Info<< "Reading gamma_a\n" << endl;

    dimensionedScalar gamma_a
    (
        transportProperties.lookup("gamma_a")
    );

Info<< "Reading gamma_h\n" << endl;

    dimensionedScalar gamma_h
    (
        transportProperties.lookup("gamma_h")
    );



 Info<< "Reading Da\n" << endl;

    dimensionedScalar Da
    (
        transportProperties.lookup("Da")
    );

 Info<< "Reading Dh\n" << endl;

    dimensionedScalar Dh
    (
        transportProperties.lookup("Dh")
    );



 Info<< "Reading aT\n" << endl;

    dimensionedScalar aT
    (
        transportProperties.lookup("aT")
    );
 
 
   
    Info<< "Reading oS\n" << endl;

    dimensionedScalar oS
    (
        transportProperties.lookup("oS")
    );

 Info<< "Reading m\n" << endl;

    dimensionedScalar m
    (
        transportProperties.lookup("m")
    );



Info<< "Reading eta\n" << endl;

    dimensionedScalar eta
    (
        transportProperties.lookup("eta")
    );



    Info<< "Reading Eb\n" << endl;

    dimensionedScalar Eb
    (
        transportProperties.lookup("Eb")
    );

    Info<< "Reading Em\n" << endl;

    dimensionedScalar Em
    (
        transportProperties.lookup("Em")
    );

    Info<< "Reading nu_b\n" << endl;

    dimensionedScalar nu_b
    (
        transportProperties.lookup("nu_b")
    );

    Info<< "Reading nu_m\n" << endl;

    dimensionedScalar nu_m
    (
        transportProperties.lookup("nu_m")
    );


    Info<< "Reading rho_b\n" << endl;

    dimensionedScalar rho_b
    (
        transportProperties.lookup("rho_b")
    );


    Info<< "Reading rho_m\n" << endl;

    dimensionedScalar rho_m
    (
        transportProperties.lookup("rho_m")
    );


    Info<< "Reading ET\n" << endl;

    dimensionedScalar ET
    (
        transportProperties.lookup("ET")
    );




    Info<< "Reading RE0\n" << endl;

    dimensionedScalar RE0
    (
        transportProperties.lookup("RE0")
    );

    Info<< "Reading b_s\n" << endl;

    dimensionedScalar b_s
    (
        transportProperties.lookup("b_s")
    );


    Info<< "Reading etaap\n" << endl;

    dimensionedScalar etaap
    (
        transportProperties.lookup("etaap")
    );


    Info<< "Reading etaapb\n" << endl;

    dimensionedScalar etaapb
    (
        transportProperties.lookup("etaapb")
    );

    Info<< "Reading tau_1\n" << endl;

    dimensionedScalar tau_1
    (
        transportProperties.lookup("tau_1")
    );


    Info<< "Reading eta_ocy\n" << endl;

    dimensionedScalar eta_ocy
    (
        transportProperties.lookup("eta_ocy")
    );

    Info<< "Reading D_T\n" << endl;

    dimensionedScalar D_T
    (
        transportProperties.lookup("D_T")
    );

    Info<< "Reading eta_apocy\n" << endl;

    dimensionedScalar eta_apocy
    (
        transportProperties.lookup("eta_apocy")
    );

    Info<< "Reading tau_2\n" << endl;

    dimensionedScalar tau_2
    (
        transportProperties.lookup("tau_2")
    );

    Info<< "Reading eta_b\n" << endl;

    dimensionedScalar eta_b
    (
        transportProperties.lookup("eta_b")
    );



        volScalarField DEV
      (
       IOobject
       (
        "DEV",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
        ),
       mesh,
       dimensionedScalar("zero",dimensionSet(0, 0, 0, 0, 0, 0, 0),0.0)
              );
	      
	volScalarField REV
      (
       IOobject
       (
        "REV",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
        ),
       DEV/(runTime.time().deltaT())
       );	      
       
        volScalarField EV
      (
       IOobject
       (
        "EV",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
        ),
       mesh,
       dimensionedScalar("zero",dimensionSet(0, 0, 0, 0, 0, 0, 0),0.0)
              );
