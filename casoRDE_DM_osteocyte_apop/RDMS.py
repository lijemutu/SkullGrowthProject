import os
import sys
import re


if not len(sys.argv)>1: #If no arguments given...
    raise Exception('Not enough arguments, given number of species and summands \n \
    eg. RDMS.py 1 4 \n\t RDMS.py 3 1 5 2')

if int(sys.argv[1])+1 != len(sys.argv[1:]): # sys.argv counts file name >python RDMS.py 1 2 3 
                                        # > ['RDMS.py', '1', '2', '3'] +2 counts N and file name
    raise Exception('Not enough summands for given number of species \n \
        eg. RDMS.py 1 4 \n\t RDMS.py 3 1 5 2')

N = int(sys.argv[1])
numSumands = [int(iCounter) for iCounter in sys.argv[2:]]

#for k in range(1,N+1): #Difussion
#    print("D_{numSpecies}" .format(numSpecies=k))

#for k in range(1,N+1): #Coefficients
#    for j in range(1,numSumands[k-1]+1):
#        print("gamma_{numSpecies}{summand}" .format(numSpecies=k,summand=j))

#for k in range(1,N+1): # Exponents
#    for j in range(1,numSumands[k-1]+1):
#        for i in range(1,N+1):
#            print("delta_{numSpecies}{summand}{product}" .format(numSpecies=k,summand=j,product=i))

class RDMSconstructor:

    def __init__(self,numberSpecies,numSumandsPerSpecie):
        """
        Initialize variables
        """
        self.N = numberSpecies
        self.numSumands = numSumandsPerSpecie
    def createFields(self):
        """
        Function to append all variables and coefficients
        needed to createFields.H of structure:

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

        and

        Info<< "Reading alpha_a\n" << endl;

        dimensionedScalar alpha_a
        (
            transportProperties.lookup("alpha_a")
        );
        """
        with open("createFields.H","a") as createFieldsFile:

            for k in range(1,self.N+1): 
                # Species definition
                speciesDefinition = "\nInfo<< \"Reading field a_{numSpecies}\\n\" << endl; \n\
                                            volScalarField a_{numSpecies}\n\
                                            ( \n\
                                                IOobject\n\
                                                    ( \n\
                                                        \"a_{numSpecies}\", \n\
                                                        runTime.timeName(), \n\
                                                        mesh, \n\
                                                        IOobject::MUST_READ, \n\
                                                        IOobject::AUTO_WRITE \n\
                                                    ), \n\
                                                mesh \n\
                                            );\n" .format(numSpecies=k)
                
                re.sub(r'(^[ \t]+|[ \t]+(?=:))', '', speciesDefinition, flags=re.M)
                createFieldsFile.write(speciesDefinition)

                speciesLocalDiffusion = "\nInfo<< \"Reading field D_{numSpecies}l\\n\" << endl; \n\
                            volScalarField D_{numSpecies}l\n\
                            ( \n\
                                IOobject\n\
                                    ( \n\
                                        \"D_{numSpecies}l\", \n\
                                        runTime.timeName(), \n\
                                        mesh, \n\
                                        IOobject::MUST_READ, \n\
                                        IOobject::AUTO_WRITE \n\
                                    ), \n\
                                mesh \n\
                            );\n" .format(numSpecies=k)

                re.sub(r'(^[ \t]+|[ \t]+(?=:))', '', speciesDefinition, flags=re.M)
                createFieldsFile.write(speciesLocalDiffusion)
            print("Species definition on createFields.H \t DONE")

            for k in range(1,self.N+1): 
                #Difussion definition
                diffDefinition = "Info<< \"Reading D_{numSpecies}\\n\" << endl; \n\
                                            dimensionedScalar D_{numSpecies} \n\
                                            ( \n\
                                                transportProperties.lookup(\"D_{numSpecies}\")\n\
                                            );\n" .format(numSpecies=k)

                re.sub(r'(^[ \t]+|[ \t]+(?=:))', '', diffDefinition, flags=re.M) 
                createFieldsFile.write(diffDefinition)
            print("Diffusion coefficients on createFields.H \t DONE")

            
            for k in range(1,self.N+1): 
                for j in range(1,self.numSumands[k-1]+1):
                    #Coefficients
                    coefficientDefinition = "Info<< \"Reading gamma_{numSpecies}{summand}\\n\" << endl; \n\
                                            dimensionedScalar gamma_{numSpecies}{summand}\n\
                                            ( \n\
                                                transportProperties.lookup(\"gamma_{numSpecies}{summand}\")\n\
                                            );\n" .format(numSpecies=k,summand=j)
                
                    re.sub(r'(^[ \t]+|[ \t]+(?=:))', '', coefficientDefinition, flags=re.M)
                   
                    createFieldsFile.write(coefficientDefinition)
            print("Coefficients definition on createFields.H \t DONE")

                    #print("gamma_{numSpecies}{summand}" .format(numSpecies=k,summand=j))

            for k in range(1,self.N+1):
                for j in range(1,self.numSumands[k-1]+1):
                    for i in range(1,self.N+1):
                        #Exponents
                        exponentDefinition = "Info<< \"Reading delta_{numSpecies}{summand}{product}\\n\" << endl; \n\
                                            dimensionedScalar delta_{numSpecies}{summand}{product}\n\
                                            ( \n\
                                                transportProperties.lookup(\"delta_{numSpecies}{summand}{product}\")\n\
                                            );\n" .format(numSpecies=k,summand=j,product=i)
                        
                        re.sub(r'(^[ \t]+|[ \t]+(?=:))', '', exponentDefinition, flags=re.M)
                        createFieldsFile.write(exponentDefinition)
            print("Expoonent coefficients on createFields.H \t DONE")

                        #print("delta_{numSpecies}{summand}{product}" .format(numSpecies=k,summand=j,product=i))
    def initialConditions(self,initialValues):
        """
        For a given number of species give a list of initial values
        2 species ---- initialConditions([0.1 , 2])
        """
        if len(initialValues)!= self.N:
            raise Exception('Not enough initial values')
        
        folderName = os.path.basename(os.path.dirname(os.path.realpath(__file__)))
        os.chdir("..")
        os.chdir("caso{folderName}" .format(folderName=folderName))
        os.chdir("0.org")

        for species in range(1,self.N+1):
            with open("a_{species}" .format(species=species) , "w") as speciesFile:
                speciesFileDefinition = "/*--------------------------------*- C++ -*----------------------------------*\\\n\
                | =========                 |                                                 |\n\
                | \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |\n\
                |  \\    /   O peration     | Version:  2.2.2                                 |\n\
                |   \\  /    A nd           | Web:      www.OpenFOAM.org                      |\n\
                |    \\/     M anipulation  |                                                 |\n\
                \\*---------------------------------------------------------------------------*/\n\
                FoamFile\n\
                {{\n\
                    version     2.0;\n\
                    format      ascii;\n\
                    class       volScalarField;\n\
                    object      a_{species};\n\
                }}\n\
                // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n\
                dimensions      [1 -3 0 0 0 0 0];\n\
                internalField	uniform {initialValue};\n" .format(species=species,initialValue=initialValues[species-1])
                speciesFile.write(speciesFileDefinition)

            with open("D_{species}l" .format(species=species), "w") as diffLocalFile:
                speciesLocalDiffusionFileDefinition = "/*--------------------------------*- C++ -*----------------------------------*\\\n\
                | =========                 |                                                 |\n\
                | \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |\n\
                |  \\    /   O peration     | Version:  2.2.2                                 |\n\
                |   \\  /    A nd           | Web:      www.OpenFOAM.org                      |\n\
                |    \\/     M anipulation  |                                                 |\n\
                \\*---------------------------------------------------------------------------*/\n\
                FoamFile\n\
                {{\n\
                    version     2.0;\n\
                    format      ascii;\n\
                    class       volScalarField;\n\
                    object      D_{species}l;\n\
                }}\n\
                // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n\
                dimensions      [1 -3 0 0 0 0 0];\n\
                internalField	uniform {initialValue};\n".format(species=species,initialValue=initialValues[species-1])

                speciesLocalDiffusionFileDefinition += 'boundaryField\n\
                {{\n\
                INNER\n\
                    {{\n\
                    type                calculated;\n\
                    value               uniform 0;\n\
                    }}\n\
                OUTER\n\
                    {{\n\
                    type                calculated;\n\
                    value               uniform 0\n\
                    }}\n\
                BOTTOM\n\
                    {{\n\
                    type                calculated;\n\
                    value               uniform 0;\n\
                    }}\n\
                WALL_INNER\n\
                    {{\n\
                    type                calculated;\n\
                    value               uniform 0;\n\
                    }}\n\
                WALL_OUTER\n\
                    {{\n\
                    type                calculated;\n\
                    value               uniform 0;\n\
                    }}\n\
                }}\n'\
                .format(species=species,initialValue=initialValues[species-1])
                diffLocalFile.write(speciesLocalDiffusionFileDefinition)

        os.chdir("..")
        os.chdir("..")
        os.chdir("{folderName}" .format(folderName=folderName))

        #print(os.listdir('.'))

    def solveDynamical(self):
        """
        Write solveDynamical.H with all summands and products
        """
        #if "solveDynamical.H" in os.listdir():
        #    os.remove("solveDynamical.H")
        
        with open("solveDynamical.H","w") as solveDynamicalFile:
            speciesHolder = ""
            
            for species in range(1,self.N+1):
                speciesHolder += "solve (\nfvm::ddt(a_{species})\n-a_{species}*fvc::div(phi)\n" .format(species=species)
                
                for summand in range(1,self.numSumands[species-1]+1):
                    speciesHolder += "-gamma_{species}{summand}*" .format(species=species,summand=summand)
                    
                    for speciesInSummand in range(1,self.numSumands[species-1]):
                        speciesHolder += "pow(a_{speciesInSummand},delta_{species}{summand}{speciesInSummand})*" \
                            .format(speciesInSummand=speciesInSummand,species=species,summand=summand)
                    speciesHolder += "pow(a_{speciesInSummand},delta_{species}{summand}{speciesInSummand})\n" \
                            .format(speciesInSummand=self.numSumands[species-1],species=species,summand=summand)
                speciesHolder+= "-fvm::laplacian(D_{species}l,a_{species}));\n\n" .format(species=species)
            solveDynamicalFile.write(speciesHolder)
        #print(speciesHolder)

    def transportProperties(self):
        """
        Function to append constants related to Diffusion coeficients, constants and exponents initialized to 0
        """
        folderName = os.path.basename(os.path.dirname(os.path.realpath(__file__)))
        os.chdir("..")
        os.chdir("caso{folderName}" .format(folderName=folderName))
        os.chdir("constant")

        with open("transportProperties.ori","a") as trasportFile:
            for k in range(1,self.N+1): #Difussion
                trasportFile.write("\nD_{numSpecies}               D_{numSpecies} [0 2 -1 0 0 0 0] 0;\n\n" .format(numSpecies=k))

            for k in range(1,self.N+1): #Coefficients
                for j in range(1,self.numSumands[k-1]+1):
                    trasportFile.write("gamma_{numSpecies}{summand}               gamma_{numSpecies}{summand} [0 0 0 0 0 0 0] 0;\n\n" .format(numSpecies=k,summand=j))

            for k in range(1,N+1): # Exponents
                for j in range(1,numSumands[k-1]+1):
                    for i in range(1,N+1):
                        trasportFile.write("delta_{numSpecies}{summand}{product}               delta_{numSpecies}{summand}{product} [0 0 0 0 0 0 0] 0;\n\n" .format(numSpecies=k,summand=j,product=i))

    def densityPoisson(self,specieDifferentiator):
        """
        Function to specify the specie needed for differentiation on mechanical properties:
        Eg: rho = (rho_m-rho_b)*(pow(oS,m)/(pow(o,m)+pow(oS,m))) + rho_b;
            rho_v = rho;
            E_v = (Em-Eb)*(pow(oS,m)/(pow(o,m)+pow(oS,m))) + Eb;
            nu_v = (nu_m-nu_b)*(pow(oS,m)/(pow(o,m)+pow(oS,m))) + nu_b;
            In this example o is the species differentiator and oS is the threshold
        """
        with open("densityPoisson.H","w") as denPoiFile:
            denPoiFile.write("rho = (rho_m-rho_b)*(pow(oS,m)/(pow(a_{specie},m)+pow(oS,m))) + rho_b;\n\
            rho_v = rho;\n\
            E_v = (Em-Eb)*(pow(oS,m)/(pow(a_{specie},m)+pow(oS,m))) + Eb;\n\
            nu_v = (nu_m-nu_b)*(pow(oS,m)/(pow(a_{specie},m)+pow(oS,m))) + nu_b;\n"\
                            .format(specie=specieDifferentiator))

    def localDiffusion(self,specieDifferentiator):
        """
        Define local diffusion properties
        Eg. 	Da_v = (pow(oS,m)/(pow(o,m)+pow(oS,m)))*Da;
	            Dh_v = (pow(oS,m)/(pow(o,m)+pow(oS,m)))*Dh;
        """
        with open("localDiffusion.H","w") as locDiffFile:
            for species in range(1,self.N+1):
                locDiffFile.write("D_{species}l =  (pow(oS,m)/(pow(a_{specieDiff},m)+pow(oS,m)))*D_{species};\n"\
                    .format(species=species,specieDiff=specieDifferentiator))

    def refineMeshField(self,specieRefinator):
        """
        Function to apply refine mesh to a given field
        Eg. oGrad = mag(fvc::grad(o));
        grad(o) is the field to refine
        """
        with open("refineMeshGuidedField.H","w") as refineFile:
            refineFile.write("oGrad = mag(fvc::grad(a_{specie}));\n" .format(specie=specieRefinator))

lee = RDMSconstructor(N,numSumands)
#lee.createFields()
#lee.initialConditions([0.1,0.2,0.3])
#lee.solveDynamical()
#lee.transportProperties()
#lee.densityPoisson(3)
#lee.localDiffusion(3)
lee.refineMeshField(3)
