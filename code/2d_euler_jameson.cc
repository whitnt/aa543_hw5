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


void readGrid(std::vector<std::vector<double> > &grid_x,
                std::vector<std::vector<double> > &grid_y,
                const int N_row, const int N_col)
{
    // Reads grid points from grid_x.dat and grid_y.dat file
    std::ifstream x_file("grid_x.dat");
    std::ifstream y_file("grid_y.dat");
    
    float x, y;
    for (int j=0; j<N_col;j++) {
        for (int i=0; i<N_row;i++) {
            x_file >> x;
            y_file >> y;
            grid_x[i][j] = x;
            grid_y[i][j] = y;
        }
    }
}

void createCellCentPoints(  std::vector<std::vector<double> > &x,
                            std::vector<std::vector<double> > &y,
                            const std::vector<std::vector<double> > grid_x,
                            const std::vector<std::vector<double> > grid_y,
                            const int N_row, const int N_col)
{
    // Uses grid points to calculate location of points at cell centers
    for (int i=0; i<N_row-1;i++) {
         for (int j=0; j<N_col-1;j++) {
            x[i][j] = 0.25*(grid_x[i][j] + grid_x[i+1][j] + grid_x[i][j+1] + grid_x[i+1][j+1]);
            y[i][j] = 0.25*(grid_y[i][j] + grid_y[i+1][j] + grid_y[i][j+1] + grid_y[i+1][j+1]);
        }
    }

    
    // add ghost cells
}

void writeCellCentPoints(const std::vector<std::vector<double> > &x,
                        const std::vector<std::vector<double> > &y,
                        const int N_row, const int N_col)
{
    // Write out coordinates of cell centered points
    std::ofstream x_file("cell_x.dat");
    std::ofstream y_file("cell_y.dat");
    
    for (int i=0; i<N_row-1;i++) {
         for (int j=0; j<N_col-1;j++) {
            x_file << x[i][j] << std::endl;
            y_file << y[i][j] << std::endl;
        }
    }
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
    // Run parameters
    
    // Free stream values
    
    // Create vectors for grid points, cell centered points, cell 
    int N_row = 129;
    int N_col = 65;
    std::vector<std::vector<double> > grid_x;
    std::vector<std::vector<double> > grid_y;
    std::vector<std::vector<double> > x;
    std::vector<std::vector<double> > y;
    std::vector<std::vector<double> > omega;
    std::vector<std::vector<double> > ds_x;
    std::vector<std::vector<double> > ds_y;
    std::vector<double> col (N_col);
    std::vector<double> colm (N_col-1);
    for (int row=0; row<N_row-1; row++) {
        grid_x.push_back(col);
        grid_y.push_back(col);
        x.push_back(colm);
        y.push_back(colm);
        omega.push_back(colm);
        ds_x.push_back(col);
        ds_y.push_back(col);
    }
    // Because the following are values at cell interfaces, they have one more 
    // Row and column than the cell valued vectors. Indices match up for the 
    // cell and it's right and lower boundary.
    grid_x.push_back(col);
    grid_y.push_back(col);
    ds_x.push_back(col);
    ds_y.push_back(col);
    
    // Read in grid and create vectors with cell centered point values,
    // Cell areas, and cell wall normals
    readGrid(grid_x, grid_y, N_row, N_col);
    createCellCentPoints(x, y, grid_x, grid_y, N_row, N_col);
    computeCellAreas();
    computeCellNormals();
    
    writeCellCentPoints(x, y, N_row, N_col);
    // Create 
    
    return 0;
}
