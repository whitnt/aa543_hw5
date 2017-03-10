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
                std::vector<std::vector<double> > &grid_y)
{
    // Reads grid points from grid_x.dat and grid_y.dat file
    int I_max = grid_x.size();
    int J_max = grid_x[0].size();
    
    std::ifstream x_file("grid_x.dat");
    std::ifstream y_file("grid_y.dat");
    
    float x, y;
    for (int j=0; j<J_max;j++) {
        for (int i=0; i<I_max;i++) {
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
                            const std::vector<std::vector<double> > grid_y)
{
    // Uses grid points to calculate location of points at cell centers
    int I_max = x.size();
    int J_max = x[0].size();
    
    for (int i=0; i<I_max;i++) {
         for (int j=0; j<J_max;j++) {
            x[i][j] = 0.25*(grid_x[i][j] + grid_x[i+1][j] + grid_x[i][j+1] + grid_x[i+1][j+1]);
            y[i][j] = 0.25*(grid_y[i][j] + grid_y[i+1][j] + grid_y[i][j+1] + grid_y[i+1][j+1]);
        }
    }

    
    // add ghost cells
}


void writeCellCentPoints(const std::vector<std::vector<double> > &x,
                        const std::vector<std::vector<double> > &y)
{
    // Write out coordinates of cell centered points
    int I_max = x.size();
    int J_max = x[0].size();
    
    std::ofstream x_file("output/cell_x.dat");
    std::ofstream y_file("output/cell_y.dat");
    
    for (int i=0; i<I_max;i++) {
         for (int j=0; j<J_max;j++) {
            x_file << x[i][j] << std::endl;
            y_file << y[i][j] << std::endl;
        }
    }
}

void computeCellAreas(std::vector<std::vector<double> > &omega,
                        const std::vector<std::vector<double> > grid_x,
                        const std::vector<std::vector<double> > grid_y
                        )
{
    // Uses grid points to calculate cell areas by takeing cross product of diags
    // Indexing of areas matches cell center points
    int I_max = omega.size();
    int J_max = omega[0].size();
    
    for (int i=0; i<I_max;i++) {
        for (int j=0; j<J_max;j++) {
            double a = grid_x[i+1][j+1] - grid_x[i][j];
            double b = grid_x[i+1][j] - grid_x[i][j+1];
            double c = grid_y[i+1][j+1] - grid_y[i][j];
            double d = grid_y[i+1][j] - grid_y[i][j+1];
            omega[i][j] = 0.5*std::abs(a*d - b*c);
        }
    }
}

void computeCellNormals(std::vector<std::vector<double> > &dsi_x,
                        std::vector<std::vector<double> > &dsi_y,
                        std::vector<std::vector<double> > &dsj_x,
                        std::vector<std::vector<double> > &dsj_y,
                        const std::vector<std::vector<double> > grid_x,
                        const std::vector<std::vector<double> > grid_y)
{
    // Uses grid points to calculate cell face normals
    // Normals are global!!
    //      dsi are normals on walls that connect verticies of increasing i
    //          and correspond to delta s_(i-1/2,j).
    //      dsj are normals on walls that connect verticies of increasing j
    //          and correspond to delta s_(i,j-1/2).
    
    int I_max = dsi_x.size();
    int J_max = dsi_x[0].size();
    
    for (int i=0; i<I_max;i++) {
        for (int j=0; j<J_max;j++) {
            dsi_x[i][j] = grid_y[i+1][j] - grid_y[i][j];
            dsi_y[i][j] = grid_x[i+1][j] - grid_x[i][j];
            dsj_x[i][j] = grid_y[i][j+1] - grid_y[i][j];
            dsj_y[i][j] = grid_x[i][j+1] - grid_x[i][j];
        }
    }
}

void setIC(std::vector<std::vector<std::vector<double> > > &u,
            const double rho_0, const double u_0, const double v_0, 
            const double E_0)
{
    // Set initial condition: u[i][j][var] = u @ infinity
    int I_max = u[0].size();
    int J_max = u[0][0].size();
    
    //~ std::cout << "I'm in the IC function!" << std::endl;
    //~ std::cout << rho_0 << u_0 << E_0 << std::endl;
    for (int i=0; i<I_max;i++) {
        for (int j=0; j<J_max;j++) {
            u[0][i][j] = rho_0;
            u[1][i][j] = rho_0*u_0;
            u[2][i][j] = rho_0*v_0;
            u[3][i][j] = rho_0*E_0;
            //~ std::cout << u[0][i][j] << u[1][i][j] << u[3][i][j] << std::endl;
        }
    }
    //~ std::cout << "number of vars = " << u[0][0].size() << std::endl;
}

void writeOutput(const std::vector<std::vector<std::vector<double> > > &u)
{
    // Write out all variables
    int I_max = u[0].size();
    int J_max = u[0][0].size();

    std::ofstream rho   ("output/rho.dat");
    std::ofstream rho_u ("output/rho_u.dat");
    std::ofstream rho_v ("output/rho_v.dat");
    std::ofstream rho_E ("output/rho_E.dat");
    
    for (int i=0; i<I_max;i++) {
        for (int j=0; j<J_max;j++) {
            rho     << u[0][i][j] << std::endl;
            rho_u   << u[1][i][j] << std::endl;
            rho_v   << u[2][i][j] << std::endl;
            rho_E   << u[3][i][j] << std::endl;
        }
    }
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
    double cfl      = 2.8;
    double R_min    = 0.000001; // Target residual value
    double vars     = 4; // Number of conserved variables
    
    // Free stream values (used in ICs and BCs)
    double gamma    = 1.4;
    double rho_0    = 0.414;
    double u_0      = 255.;
    double v_0      = 0.;
    double c_0      = 300.;
    double E_0      = 193.;
    
    // Allocate all global vectors that will be used throughout sim
    int N_row = 129;
    int N_col = 65;
    // Grid points at cell verticies
    std::vector<std::vector<double> > grid_x;
    std::vector<std::vector<double> > grid_y;
    // Grid points at cell centers
    std::vector<std::vector<double> > x;
    std::vector<std::vector<double> > y;
    // Cell areas
    std::vector<std::vector<double> > omega;
    // Cell wall normals
    std::vector<std::vector<double> > dsi_x;
    std::vector<std::vector<double> > dsi_y;
    std::vector<std::vector<double> > dsj_x;
    std::vector<std::vector<double> > dsj_y;
    // Conserved variables
    std::vector<std::vector<std::vector<double> > > u;
    // Size the vectors
    std::vector<double> col (N_col); // For vectors that hold verticie values
    std::vector<double> colm (N_col-1); // For vectors that hold cell centered values
    grid_x.push_back(col);
    grid_y.push_back(col);
    for (int row=0; row<N_row-1; row++) {
        grid_x.push_back(col);
        grid_y.push_back(col);
        x.push_back(colm);
        y.push_back(colm);
        omega.push_back(colm);
        dsi_x.push_back(colm);
        dsi_y.push_back(colm);
        dsj_x.push_back(colm);
        dsj_y.push_back(colm);
    }
    for (int var=0;var<vars;var++){
        u.push_back(x);
    }
    
    // Prep work:
    //      Read in grid with cell verticies
    //      Find cell centers
    //      Find cell areas
    //      Find cell wall normals
    readGrid(grid_x, grid_y);
    createCellCentPoints(x, y, grid_x, grid_y);
    computeCellAreas(omega, grid_x, grid_y);
    computeCellNormals(dsi_x, dsi_y, dsj_x, dsj_y, grid_x, grid_y);
    
    writeCellCentPoints(x, y);
    
    // Set up initial conditions
    setIC(u, rho_0, u_0, v_0, E_0);
    writeOutput(u);
    
    // Run sim
    
    return 0;
}
