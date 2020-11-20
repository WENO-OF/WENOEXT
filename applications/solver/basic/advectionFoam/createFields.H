	 Info<< "Reading field T\n" << endl;
	 
	 volScalarField T
	 (
		 IOobject
		 (
			 "T",
			 runTime.timeName(),
			 mesh,
			 IOobject::MUST_READ,
			 IOobject::AUTO_WRITE
		 ),
		 mesh
	 );
	 
	 
	 Info<< "Reading field U\n" << endl;
	 
	 volVectorField U
	 (
		 IOobject
		 (
			 "U",
			 runTime.timeName(),
			 mesh,
			 IOobject::MUST_READ,
			 IOobject::AUTO_WRITE
		 ),
		 mesh
	 );
	 
	 
	 #include "createPhi.H"
 

	singlePhaseTransportModel laminarTransport(U, phi);

    autoPtr<incompressible::RASModel> turbulence
    (
        incompressible::RASModel::New(U, phi, laminarTransport)
    );
