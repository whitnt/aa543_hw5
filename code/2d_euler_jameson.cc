/*-------------------------------------------------------

* 2D Euler equations solved for steady state condition on 
* non-uniform grid using the Jameson scheme


---------------------------------------------------------*/

// Import from std
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <vector>
#include <cstdlib>


void readGrid()
{
    // Reads grid points from grid.dat file

}

void createCellCentPoints()
{
    // Uses grid points to calculate location of points at cell centers
    // x[i][j] = 0.25*(grid_x[i][j] + grid_x[i+1][j] + grid_x[i][j+1] + grid_x[i+1][j+1])
    // y[i][j] = 0.25*(grid_y[i][j] + grid_y[i+1][j] + grid_y[i][j+1] + grid_y[i+1][j+1])
    
    // add ghost cells
}

void computeCellAreas()
{
    // Uses grid points to calculate cell areas 
    // Indexing of areas matches cell center points
    // omega[i][j] = 
}

void computeCellNormals()
{
    // Uses grid points to calculate cell face normals
    //      Face 0 -> dsx[i][j][0] = 
    //      Face 1 -> dsy[i][j][1] = 
}

void setIC()
{
    // Set initial condition: u[i][j][var] = u @ infinity
    //  u[i][j][0] = rho_0;
    //  u[i][j][1] = rho_0*u_0;
    //  u[i][j][2] = rho_0*v_0;
    //  u[i][j][3] = rho_0*E_0;
}

void writeOutput()
{
    
}

void RK1()
{
    
}

void setBC()
{
    
}

int main()
{
    // 

    return 0;
}
