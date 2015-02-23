/*---------------------------------------------------------------------------*\
Date started: 10/05/2013

Moves the grid points.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "mathematicalConstants.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"

    IOdictionary initDict
    (
        IOobject
        (
            "ScharMountainDict",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );
    
    IOField<point> newPoints
    (
        IOobject("points", mesh.time().constant(), "polyMesh", mesh),
        mesh.points()
    );

    // Get which coord system to use
    const scalar coordSys(readScalar(initDict.lookup("coordSys")));
    
    // Declare constants
    const scalar zt(readScalar(initDict.lookup("zt"))); // domain height
    const scalar a(readScalar(initDict.lookup("a")));
    const scalar lam(readScalar(initDict.lookup("lam")));
    const scalar hm(readScalar(initDict.lookup("hm")));
    
    // Declare height of ground
    scalarField h(newPoints.size());
    scalarField hStar(newPoints.size());
    scalarField h1(newPoints.size());
    scalarField h2(newPoints.size());
    const scalar zh(readScalar(initDict.lookup("zh")));
    const scalar s(readScalar(initDict.lookup("s")));
    const scalar s1(readScalar(initDict.lookup("s1")));
    const scalar s2(readScalar(initDict.lookup("s2")));
    const scalar n(readScalar(initDict.lookup("n")));
    float b1,b2,beta;
    int posn(0);
    
    
    // create height field
    for(label ip = 0; ip < newPoints.size(); ip++)
    {   
        if (abs(newPoints[ip].x())<=a)
            {
            hStar[ip] = hm*pow(Foam::cos(M_PI*newPoints[ip].x()/(2*a)),2);
            
            h[ip] =pow(Foam::cos(M_PI*newPoints[ip].x()/(lam)),2)*hStar[ip];
            }
        else { h[ip]=0; }
        
        if (coordSys == 2)
        {
            // SLEVE heights
            h1[ip] = 0.5*hStar[ip]; // smoothed terrain
            h2[ip] = h[ip]-h1[ip]; // remainder
        }
    }
        
    // move all grid points
    for(label ip = 0; ip < newPoints.size(); ip++)
    {   
        if (abs(newPoints[ip].x())<=a)
        {
            h[ip] =pow(Foam::cos(M_PI*newPoints[ip].x()/(lam)), 2)
                  *hm*pow(Foam::cos(M_PI*newPoints[ip].x()/(2*a)), 2);
        }
        else { h[ip]=0; }   
     
        // move the z points
        
        if (coordSys == 0)
        {
            // BTF
            newPoints[ip].z() = newPoints[ip].z()*(zt-h[ip])/zt + h[ip];
        }
        
        if (coordSys == 1)
        {
            // HTF 
            newPoints[ip].z() = newPoints[ip].z()
                 + h[ip]*Foam::sinh((zt-newPoints[ip].z())/s)/Foam::sinh(zt/s);
        }
        
        if (coordSys == 2)
        {
            // SLEVE
            b1 = Foam::sinh
            (
                (Foam::pow(zt/s1,n)-Foam::pow(newPoints[ip].z()/s1,n))
            )/Foam::sinh(Foam::pow(zt/s1,n));
            
            b2 = Foam::sinh(Foam::pow(zt/s2,n)-Foam::pow(newPoints[ip].z()/s2,n))
                /Foam::sinh(Foam::pow(zt/s2,n));
                
            newPoints[ip].z() = newPoints[ip].z() + h1[ip]*b1+h2[ip]*b2;
        }
              
//        // will save to a file then plottable by matlab 
//        newPts << newPoints[ip].x() << " " << newPoints[ip].y() << " "
//              << newPoints[ip].z() << nl;
               
        if (coordSys == 3)
        {
	
            // !!! STF ONLY !!!! if moved up a level apply move smoothing for the height field
            // if the previous y value is 1000 and the current is 0 then 
            // the values belong to the next line up
            if (ip>2 && newPoints[ip].y()==0 && newPoints[ip-1].y()!=0)
            {
                //beta = 0.2*Foam::min(newpoints[posnip].z());
                //h[ip] = h[ip-1] + beta;
                //posn = ip;
            }
        }
    }
    
    newPoints.write();
}
