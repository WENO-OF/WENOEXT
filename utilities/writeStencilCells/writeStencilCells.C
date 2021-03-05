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
    
    argList::addOption
    (
        "cellIndex",
        "label",
        "give the cell index for the stencils to be written"
    );
    
    #include "setRootCase.H"
    #include "createTime.H"        // create the time object
    #include "createMesh.H"        // create the mesh object

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    
    
    label polOrder = 3;
    label cellIndex = 0;
    args.optionReadIfPresent("polOrder", polOrder);
    args.optionReadIfPresent("cellIndex", cellIndex);
    
    const WENOBase& WENO = WENOBase::instance(mesh,polOrder);
    
    // Get the cell stencil list
    auto& stencilID = WENO.stencilsID()[cellIndex];
    
    forAll(stencilID,stencilI)
    {
        // Create a field for invalid cells 
        volScalarField field 
        (
            IOobject
                (
                    "Cell-"+name(cellIndex)+"-stencil"+name(stencilI),
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh,
                dimensionedScalar("0",dimless,-1)
        );
        
        forAll(stencilID[stencilI],cellJ)
        {
            if (stencilID[stencilI][cellJ] != int(WENOBase::Cell::deleted)
                //|| stencilID[stencilI][cellJ] != int(WENOBase::Cell::empty)
                )
            field[stencilID[stencilI][cellJ]] = stencilI;
        }
        field[cellIndex] = 100;
        field.write();
    }

    return 0;
}

