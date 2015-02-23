/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2010-2010 OpenCFD Ltd.
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
    testEigenSVD

Description
    Tests the Eigen library SVD analysis

\*---------------------------------------------------------------------------*/

#include "IOstreams.H"
#include "OSspecific.H"
#include "OFstream.H"

#include <iostream>
#include <Eigen/SVD>

//using Eigen::MatrixXd;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

using namespace Foam;
//using namespace std;
//using namespace Eigen;

int main(int argc, char *argv[])
{
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    // The matrix to decompose
    Eigen::MatrixXd M(6,3);
    M(0,0) = 1e6;
    M(0,1) = -1e6;
    M(0,2) = 0;
    M(1,0) = 1e6;
    M(1,1) = 1e6;
    M(1,2) = 0;
    M(2,0) = 1e3;
    M(2,1) = -1e3;
    M(2,2) = 0.6;
    M(3,0) = 1e3;
    M(3,1) = 1e3;
    M(3,2) = 0.6;
    M(4,0) = 1e3;
    M(4,1) = -1e3;
    M(4,2) = -1.2;
    M(5,0) = 1e3;
    M(5,1) = 1e3;
    M(5,2) = -1.2;

    Eigen::JacobiSVD<Eigen::MatrixXd> svd(M, Eigen::ComputeThinU | Eigen::ComputeThinV);
    cout << "Singular values\n" << svd.singularValues() << '\n';
    Eigen::MatrixXd rhs = Eigen::MatrixXd::Identity(M.rows(), M.rows());
    cout << "\nrhs = \n" << rhs << '\n';
    cout << "\nweights = \n" << svd.solve(rhs) << '\n';
//cout << "Its left singular vectors are the columns of the thin U matrix:" << endl << svd.matrixU() << endl;
//cout << "Its right singular vectors are the columns of the thin V matrix:" << endl << svd.matrixV() << endl;
//Vector3f rhs(1, 0, 0);
//cout << "Now consider this rhs vector:" << endl << rhs << endl;
//cout << "A least-squares solution of m*x = rhs is:" << endl << svd.solve(rhs) << endl;    


    Info<< "\nEnd\n" << endl;
    return 0;
}


// ************************************************************************* //
