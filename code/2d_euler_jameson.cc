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
    for (int j = 0; j < J_max; j++) {
        for (int i = 0; i < I_max; i++) {
            x_file >> x;
            y_file >> y;
	    if (i == 0) {
		std::cout << "x is" << x << "\n";
		std::cout << "y is" << y << "\n";
	    }
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
    
    for (int i = 0; i < I_max; i++) {
         for (int j = 0; j < J_max; j++) {
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
    
    for (int i=0; i < I_max;i++) {
         for (int j=0; j < J_max;j++) {
            x_file << x[i][j] << std::endl;
            y_file << y[i][j] << std::endl;
        }
    }
}


void computeCellAreas(std::vector<std::vector<double> > &omega,
                        const std::vector<std::vector<double> > grid_x,
                        const std::vector<std::vector<double> > grid_y)
{
    // Uses grid points to calculate cell areas by taking cross product of diags
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
    
    for (int i = 0; i < I_max; i++) {
        for (int j = 0; j < J_max; j++) {
	  dsi_x[i][j] = (grid_y[i+1][j] - grid_y[i][j]);
	  dsi_y[i][j] = -1.0*(grid_x[i+1][j] - grid_x[i][j]); // Modified w/ -1
	  dsj_x[i][j] = (grid_y[i][j+1] - grid_y[i][j]);
	  dsj_y[i][j] = -1.0*(grid_x[i][j+1] - grid_x[i][j]); // Modified w/ -1
        }
    }
}

void setIC(std::vector< std::vector< std::vector<double> > > &u,
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

void spaceInt()
{
    // 
    
}

void calcTau(const std::vector<std::vector<std::vector<double> > > &u, double tau_ij)
{
    // Calculate local time step
    // Given parameters
    double cfl = 2.8;
    double gamma = 1.4;
    
    // Find local u, v and c
    double rho_ij = u[0][i][j];
    double u_ij = u[1][i][j]/rho_ij;
    double v_ij = u[2][i][j]/rho_ij;
    double p_ij = (gamma - 1)*(0.5*rho_ij*std::sqrt(u_ij*u_ij + v_ij*v_ij) 
                            - u[3][i][j]);
    double c_ij = std::sqrt(gamma*p_ij/rho_ij);
    
    // Average wall normals
    double dsi_x_avg = 0.5*(dsi_x[i][j] + dsi_x[i][j+1]);
    double dsi_y_avg = 0.5*(dsi_y[i][j] + dsi_y[i][j+1]);
    double dsj_x_avg = 0.5*(dsi_x[i][j] + dsi_x[i+1][j]);
    double dsj_y_avg = 0.5*(dsi_y[i][j] + dsi_y[i+1][j]);
    
    // Calculate max possible speed
    double uc = u_ij + c_ij;
    double vc = v_ij + c_ij;
    
    // Calculate time step
    tau_ij = cfl/(std::abs(dsi_x_avg*uc + dsi_y_avg*vc) 
                                + std::abs(dsj_x_avg*uc + dsj_y_avg*vc));
}

void tempInt(const std::vector<std::vector<std::vector<double> > > &u,
            const std::vector<std::vector<std::vector<double> > > &r,
            const std::vector<std::vector<double> > &dsi_x,
            const std::vector<std::vector<double> > &dsi_y,
            const std::vector<std::vector<double> > &dsj_x,
            const std::vector<std::vector<double> > &dsj_y,
            const double alpha,
            std::vector<std::vector<std::vector<double> > > &u_new)
{
    // Computes RK step given variable vector u, residual from spacial 
    // integrator r, a local time step (calculated here) and RK step weight alpha.
    
    double tau_ij = 0.;
    for (int i; i<u[0].size; i++){
        for (int j; j<u[0][0].size(); j++) {
            // Calculate local time step
            calcDT(u, tau_ij);
            
            // Calc RK step
            u_new[0][i][j] = u[0][i][j] - alpha*tau_ij*r[0][i][j];
            u_new[1][i][j] = u[1][i][j] - alpha*tau_ij*r[1][i][j];
            u_new[2][i][j] = u[2][i][j] - alpha*tau_ij*r[2][i][j];
            u_new[3][i][j] = u[3][i][j] - alpha*tau_ij*r[3][i][j];
        }
    }
}

void setBC()
{
    
}

int main()
{
    // Run parameters
    double cfl      = 2.8;
    double R_min    = 1.0e-6; // Target residual value
    double vars     = 4; // Number of conserved variables
    
    // Free stream values (used in ICs and BCs)
    double gamma    = 1.4;
    double rho_0    = 0.414;
    double u_0      = 255.;
    double v_0      = 0.;
    double c_0      = 300.;
    double E_0      = 193.;
    
    ////// Allocate all global vectors that will be used throughout sim
    int N_row = 129;
    int N_col = 65;
    // Grid points at cell verticies
    std::vector< std::vector<double> > grid_x;
    std::vector< std::vector<double> > grid_y;
    // Grid points at cell centers
    std::vector< std::vector<double> > x;
    std::vector< std::vector<double> > y;
    // Cell areas
    std::vector< std::vector<double> > omega;
    // Cell wall normals
    std::vector< std::vector<double> > dsi_x;
    std::vector< std::vector<double> > dsi_y;
    std::vector< std::vector<double> > dsj_x;
    std::vector< std::vector<double> > dsj_y;
    // Conserved variables
    std::vector< std::vector< std::vector<double> > > u;
    
    ////// Size the vectors
    std::vector<double> col (N_col); // For vectors that hold verticie values
    std::vector<double> colm (N_col-1); // For vectors that hold cell centered values
    grid_x.push_back(col);
    grid_y.push_back(col);
    
    for (int row = 0; row < N_row - 1; row++) {
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
    
    for (int var = 0; var<vars; var++){
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
    
    // Run sim (Using standard RK4 in time and Jameson artifical viscosty in space)
    double R = 1.0;
    while (R > 1.0e-6) {
        // Create copies of u for RK steps, plus an r to store residuals;
        std::vector< std::vector< std::vector<double> > > u_1(u);
        std::vector< std::vector< std::vector<double> > > u_2(u);
        std::vector< std::vector< std::vector<double> > > u_3(u);
        std::vector< std::vector< std::vector<double> > > u_4(u);
        std::vector< std::vector< std::vector<double> > > r(u);
        
        // First RK step
        // Calculate residual
        spaceInt(u, omega, dxi_x, dsi_y, dsj_x, dsj_y, r); 
        // Calculate first RK step
        double alpha = 0.25;
        tempInt(u, r, omega, dxi_x, dsi_y, dsj_x, dsj_y, alpha, u_1);
        
        // Second RK step
        // Calculate residual using u_1
        spaceInt(u_1, omega, dxi_x, dsi_y, dsj_x, dsj_y, r); 
        // Calculate RK step
        alpha = 1./3.;
        tempInt(u, r, omega, dxi_x, dsi_y, dsj_x, dsj_y, alpha, u_2);
        
        // Third RK step
        // Calculate residual using u_2
        spaceInt(u_2, omega, dxi_x, dsi_y, dsj_x, dsj_y, r); 
        // Calculate RK step
        alpha = 0.5;
        tempInt(u, r, omega, dxi_x, dsi_y, dsj_x, dsj_y, alpha, u_3);
        
        // Fourth RK step
        // Calculate residual using u_3
        spaceInt(u_3, omega, dxi_x, dsi_y, dsj_x, dsj_y, r); 
        // Calculate RK step
        alpha = 1.;
        tempInt(u, r, omega, dxi_x, dsi_y, dsj_x, dsj_y, alpha, u_4);
        
        // Apply boundary conditions
        applyBC(u_4);
        
        // Update original vector u
        
        // Calculate total residual
        
        
    }
    
    writeOutput(u);
    
    return 0;
}
