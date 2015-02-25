#include "fvCFD.H"

class VelocityProfile
{
    public:
    VelocityProfile(const IOdictionary& dict)
        :
            u0(readScalar(dict.lookup("maxVelocity"))),
            z1(readScalar(dict.lookup("zeroVelocityHeight"))),
            z2(readScalar(dict.lookup("maxVelocityHeight")))
    {}

    template<class Type, template<class> class PatchField, class GeoMesh>
    void applyTo(GeometricField<Type, PatchField, GeoMesh>& field)
    {
        applyToInternalField(field);
        applyToBoundary("inlet", field);
        applyToBoundary("outlet", field);
    }

    private:
    const scalar u0;
    const scalar z1;
    const scalar z2;
    
    template<class Type>
    void applyToInternalField(GeometricField<Type, fvsPatchField, surfaceMesh>& field)
    {
        forAll(field, faceI)
        {
            const point& Cf = field.mesh().Cf()[faceI];
            field[faceI] = velocityAt(Cf.z());
        }
    }

    template<class Type>
    void applyToInternalField(GeometricField<Type, fvPatchField, volMesh>& field)
    {
        forAll(field, cellI)
        {
            const point& C = field.mesh().C()[cellI];
            field[cellI] = velocityAt(C.z());
        }
    }

    template<class Type, template<class> class PatchField, class GeoMesh>
    void applyToBoundary
    (
        const word name, GeometricField<Type, PatchField, GeoMesh>& field
    )
    {
        label boundaryI = findBoundaryPatchIndex(field.mesh(), name);
        forAll(field.boundaryField()[boundaryI], cellI)
        {
            const point& face = field.mesh().Cf().boundaryField()[boundaryI][cellI];
            field.boundaryField()[boundaryI][cellI] = velocityAt(face.z());
        }
    }

    vector velocityAt(const scalar z) const
    {
        if (z > z1 && z < z2)
        {
            return vector(u0*sqr(Foam::sin(M_PI/2*(z-z1)/(z2-z1))), 0, 0);
        }
        else if (z >= z2)
        {
            return vector(u0, 0, 0);
        }
        else
        {
            return vector(0, 0, 0);
        }
    }

    label findBoundaryPatchIndex(const fvMesh& mesh, const word& name)
    {
        forAll(mesh.boundaryMesh(), patchI)
        {
            if (mesh.boundaryMesh()[patchI].name() == name)
            {
                return patchI;
            }
        }
        FatalErrorIn("setHorizontalVelocityField")
            << " no boundary called " << name << ". The boundaries are called "
            << mesh.boundaryMesh().names()
            << exit(FatalError);

        return -1;
    }
};

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"
    Info << "Reading velocityFieldDict" << endl;

    IOdictionary initDict
    (
        IOobject
        (
            "velocityFieldDict",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    surfaceVectorField Uf
    (
        IOobject("Uf", runTime.timeName(), mesh),
        mesh,
        dimensionedVector("Uf", dimVelocity, vector(0,0,0)),
        "fixedValue"
    );

    VelocityProfile velocityProfile(initDict);
    Info << "Creating velocity field Uf" << endl;
    velocityProfile.applyTo(Uf);
    Uf.write();

    volVectorField U
    (
        IOobject("U", runTime.timeName(), mesh),
        mesh,
        dimensionedVector("U", dimVelocity, vector(0,0,0)),
        "fixedValue"
    );

    Info << "Creating velocity field U" << endl;
    velocityProfile.applyTo(U);

    U.write();
}
