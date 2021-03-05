#include "fvCFD.H"                 // include basic openFoam classes
#include "WENOBase.H"

int main(int argc, char *argv[])   // start main loop
{
    argList::addOption
    (
        "polOrder",
        "label",
        "Specify the polOrder of the WENO scheme"
        "(default 3)"
    );
    #include "setRootCase.H"
    #include "createTime.H"        // create the time object
    #include "createMesh.H"        // create the mesh object

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    
    
    label polOrder = 3;
    args.optionReadIfPresent("polOrder", polOrder);
    
    const WENOBase& WENO = WENOBase::instance(mesh,polOrder);
    
    // Create a field for invalid cells 
    volScalarField field 
    (
        IOobject
            (
                "deletedCells.WENO",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("0",dimless,0)
    );
    
    label invalidCells = 0;

    forAll(WENO.stencilsID(),cellI)
    {
        if (WENO.stencilsID()[cellI][0][0] == int(WENOBase::Cell::empty))
        {
            field[cellI] = 1;
            invalidCells++;
        }
    }

    field.write();
    Info << "Number of deleted cells: "<<invalidCells<<endl;
    
    // ------------------------------------------------------------------------
    //                  Percentage of valid stencils per cell
    // ------------------------------------------------------------------------
    
    // Create a field for invalid cells 
    volScalarField validStencils 
    (
        IOobject
            (
                "validStencilRate.WENO",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("0",dimless,0)
    );
    

    forAll(WENO.stencilsID(),cellI)
    {
        const label numStencil = WENO.stencilsID()[cellI].size();
        label invalidStencils = 0;
        forAll(WENO.stencilsID()[cellI],stencilI)
        {
            if (WENO.stencilsID()[cellI][stencilI][0] == int(WENOBase::Cell::deleted))
                invalidStencils++;
        }
        validStencils[cellI] = 1.0 - double(invalidStencils)/double(numStencil);
    }
    validStencils.write();
    
    
    return 0;
}

