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
#include <algorithm>
#include <vector>
#include <cstdlib>


void readGrid(std::vector<std::vector<double> > &grid_x,
                std::vector<std::vector<double> > &grid_y)
{
    // Reads grid points from grid_x.dat and grid_y.dat file
    int I_max = grid_x.size();
    int J_max = grid_x[0].size();
    
    std::ifstream x_file("output/grid_x.dat");
    std::ifstream y_file("output/grid_y.dat");
    
    float x, y;
    for (int j = 0; j < J_max; j++) {
        for (int i = 0; i < I_max; i++) {
            x_file >> x;
            y_file >> y;
	    // For debug
	    //if (i == 0) {
	    //std::cout << "x is" << x << "\n";
	    //	std::cout << "y is" << y << "\n";
	    //}
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
    
    for (int i = 0; i < x.size(); i++) {
         for (int j = 0; j < x[0].size(); j++) {
            x[i][j] = 0.25*(grid_x[i][j] + grid_x[i+1][j] + grid_x[i][j+1] + grid_x[i+1][j+1]);
            y[i][j] = 0.25*(grid_y[i][j] + grid_y[i+1][j] + grid_y[i][j+1] + grid_y[i+1][j+1]);
        }
    }
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
    
    // Do dsi_x and dsi_y first
    for (int i = 0; i < dsi_x.size(); i++) {
        for (int j = 0; j < dsi_x[0].size(); j++) {
        // Calculate them
        dsi_x[i][j] = -1.0*(grid_y[i+1][j] - grid_y[i][j]);
        dsi_y[i][j] = 1.*(grid_x[i+1][j] - grid_x[i][j]); // Modified w/ -1
        }
    }
    
    // then do dsj_x and dsj_y
    for (int i = 0; i < dsj_x.size(); i++) {
        for (int j = 0; j < dsj_x[0].size(); j++) {
        // Calculate them
        dsj_x[i][j] = (grid_y[i][j+1] - grid_y[i][j]);
        dsj_y[i][j] = -1.0*(grid_x[i][j+1] - grid_x[i][j]); // Modified w/ -1
        }
    }
}

void writeCellNormals(const std::vector<std::vector<double> > &dsi_x,
                      const std::vector<std::vector<double> > &dsi_y,
                      const std::vector<std::vector<double> > &dsj_x,
                      const std::vector<std::vector<double> > &dsj_y)
{
    // Write out coordinates of cell centered points
    std::ofstream ix_file("output/dsi_x.dat");
    std::ofstream iy_file("output/dsi_y.dat");
    std::ofstream jx_file("output/dsj_x.dat");
    std::ofstream jy_file("output/dsj_y.dat");
    
    for (int i=0; i < dsi_x.size();i++) {
         for (int j=0; j < dsi_x[0].size();j++) {
            ix_file << dsi_x[i][j] << std::endl;
            iy_file << dsi_y[i][j] << std::endl;
        }
    }

    for (int i=0; i < dsj_x.size();i++) {
         for (int j=0; j < dsj_x[0].size();j++) {
            jx_file << dsj_x[i][j] << std::endl;
            jy_file << dsj_y[i][j] << std::endl;
        }
    }
}

void setIC(std::vector< std::vector< std::vector<double> > > &u,
            const double rho_0, const double u_0, const double v_0, 
            const double p_0, const double gamma)
{
    // Set initial condition: u[i][j][var] = u @ infinity
    
    //~ std::cout << "I'm in the IC function!" << std::endl;
    //~ std::cout << rho_0 << u_0 << E_0 << std::endl;
    for (int i=0; i<u[0].size();i++) {
        for (int j=0; j<u[0][0].size();j++) {
            u[0][i][j] = rho_0;
            u[1][i][j] = rho_0*u_0;
            u[2][i][j] = rho_0*v_0;
            u[3][i][j] = 1./(gamma - 1.)*p_0 + rho_0*(u_0*u_0 + v_0*v_0)/2.;
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
    std::ofstream ghost_ext ("output/ghost_ext.dat");
    std::ofstream ghost_int ("output/ghost_air.dat");
    
    for (int i=0; i<I_max;i++) {
        for (int j=0; j<J_max-2;j++) {
            rho     << u[0][i][j] << std::endl;
            rho_u   << u[1][i][j] << std::endl;
            rho_v   << u[2][i][j] << std::endl;
            rho_E   << u[3][i][j] << std::endl;
            
        }
    }
                
    for (int i=0; i<I_max;i++) { 
        ghost_ext << u[0][i][u[0][0].size() - 1] << "\t\t" 
                << u[1][i][u[0][0].size() - 1] << "\t\t" 
                << u[2][i][u[0][0].size() - 1] << "\t\t"
                << u[3][i][u[0][0].size() - 1] << "\t\t"
                << std::endl;
    }
    
    for (int i=0; i<I_max;i++) { 
        ghost_int << u[0][i][u[0][0].size() - 2] << "\t\t" 
                << u[1][i][u[0][0].size() - 2] << "\t\t" 
                << u[2][i][u[0][0].size() - 2] << "\t\t"
                << u[3][i][u[0][0].size() - 2] << "\t\t"
                << std::endl;
    }
}

void spaceInt(const std::vector<std::vector<std::vector<double> > > &u,
            const std::vector<std::vector<double> > &omega,
            const std::vector<std::vector<double> > &dsi_x,
            const std::vector<std::vector<double> > &dsi_y,
            const std::vector<std::vector<double> > &dsj_x,
            const std::vector<std::vector<double> > &dsj_y,
	    std::vector<std::vector<std::vector<double> > > &r,
	    const double &gamma)
{
    // Calculate residual r using Jameson scheme with artificial viscosity
    // Initialize stuff
    std::vector<std::vector<std::vector<double> > > F(u);
    std::vector<std::vector<std::vector<double> > > G(u);
    std::vector<std::vector<double> > p(omega);
    std::vector<std::vector<double> > sound(omega);
    std::vector<std::vector<double> > nui(omega);
    std::vector<std::vector<double> > nuj(omega);

    // Initialize "D" viscosity vectors
    //~ std::vector<std::vector<std::vector<double> > > dpi12(omega);
    //~ std::vector<std::vector<std::vector<double> > > dmi12(omega);
    //~ std::vector<std::vector<std::vector<double> > > dpj12(omega);
    //~ std::vector<std::vector<std::vector<double> > > dmj12(omega);
    
    // Calculate and store cell centered fluxes
    for (int i=0; i < u[0].size(); i++){
        for (int j=0; j < u[0][0].size(); j++) {
            // Create temp variables for readability
            double u_0 = u[0][i][j];
            double u_1 = u[1][i][j];
            double u_2 = u[2][i][j];
            double u_3 = u[3][i][j];
            double p_ij = 0.;
            
            if (j >= u[0][0].size()-2) {
                p_ij = (gamma - 1.)*(u_3 - 0.5*(u_1*u_1 + u_2*u_2)/u_0);
            } else {
                p[i][j] = (gamma - 1.)*(u_3 - 0.5*(u_1*u_1 + u_2*u_2)/u_0);
                sound[i][j] = sqrt(gamma*p[i][j]/u_0);
                p_ij = p[i][j];
            }
            
            // Calculate fluxes
            F[0][i][j] = u_1;
            //~ std::cout << F[0][i][j] << std::endl;
            F[1][i][j] = u_1*u_1/u_0 + p_ij;
            F[2][i][j] = u_1*u_2/u_0;
            F[3][i][j] = u_1*(u_3 + p_ij)/u_0;
            G[0][i][j] = u_2;
            G[1][i][j] = u_1*u_2/u_0;
            G[2][i][j] = u_2*u_2/u_0 + p_ij;
            G[3][i][j] = u_2*(u_3 + p_ij)/u_0;
            
        }
    }
    
    
    //~ writeOutput(F);
    
    //~ double eps2pi12(omega);
    //~ double eps2mi12(omega);
    //~ double eps2pj12(omega);
    //~ double eps2mj12(omega);
    //~ double eps4pi12(omega);
    //~ double eps4mi12(omega);
    //~ double eps4pj12(omega);
    //~ double eps4mj12(omega);
    
    ////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////
    /////////////////////////NORMALIZED PRESSURE GRADIENT //////////////////////////
    ////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////
    // Calculate the "nu" values for epsilon viscosity parameters
    // Normalized pressure gradient
    // NEEDS: Boundary Info
    //~ for (int i = 1; i < (u[0].size() - 1); i++) {
      //~ for (int j = 1; j < (u[0][0].size() - 1); j++) {
	//~ nui[i][j] = abs(p[i+1][j] - 2.0*p[i][j] + p[i-1][j])/
	  //~ (p[i+1][j] + 2.0*p[i][j] + p[i-1][j]);
	//~ nuj[i][j] = abs(p[i][j+1] - 2.0*p[i][j] + p[i][j-1])/
	  //~ (p[i][j+1] + 2.0*p[i][j] + p[i][j-1]);
      //~ }
    //~ }
    
    //~ int iSZ = u[0].size() - 1; // For indexing in periodic/ghost fixes
    //~ int jSZ = u[0][0].size() - 3; // For indexing in periodic/ghost fixes
    //~ // Calculate the D viscosity parameters for cell interfaces
    //~ for (int i = 0; i < u[0].size(); i++) {
      //~ for (int j = 0; (j < u[0][0].size() - 2); j++) {
	//~ double u_ij = u[1][i][j]/u[0][i][j];
	//~ double v_ij = u[2][i][j]/u[0][i][j];
	//~ //////////////////////////////////////////////////////////////////////////////
	//~ ///////////////// CALCULATION FOR I'S ////////////////////////////////////////
	//~ ///////////////////////////////////////////////// CALCULATION FOR DPI12, DMI12
	//~ //////////////////////////////////////////////////////////////////////////////
	//~ //////////////////////////////////////////////////////////////////////////////
	//~ //////////////////////////////////////////////////////////////////////////////
	//~ // Calculate epsilon parameters
	//~ // Epsilon2 calc
	//~ if ( i == 0 ) { // Periodic "i" Fix
	  //~ // Wave Speed PLUS
	  //~ double wave_ip12 = u_ij*dsj_x[i+1][j] + v_ij*dsj_y[i+1][j] +
	    //~ sound[i][j]*sqrt(dsj_x[i+1][j]*dsj_x[i+1][j] + dsj_y[i+1][j]*dsj_y[i+1][j]);
	  //~ // Wave Speed MINUS
	  //~ double wave_im12 = -u_ij*dsj_x[i][j] - v_ij*dsj_y[i][j] +
	    //~ sound[i][j]*sqrt(dsj_x[i][j]*dsj_x[i][j] + dsj_y[i][j]*dsj_y[i][j]);
	  //~ // Epsilon12 PLUS
	  //~ double eps2pi12 = 0.5*0.25*wave_ip12*
	    //~ std::max(nui[iSZ][j], nui[i][j], nui[i+1][j], nui[i+2][j]);
	  //~ // Epsilon12 MINUS
	  //~ double eps2mi12 = 0.5*0.25*wave_im12*
	    //~ std::max(nui[iSZ][j], nui[i][j], nui[i+1][j], nui[i+2][j]);
	//~ }
	//~ if ( i == u[0].size() - 2 ) { // Periodic "i" Fix
	  //~ // Wave speed PLUS
	  //~ double wave_ip12 = u_ij*dsj_x[i+1][j] + v_ij*dsj_y[i+1][j] +
	    //~ sound[i][j]*sqrt(dsj_x[i+1][j]*dsj_x[i+1][j] + dsj_y[i+1][j]*dsj_y[i+1][j]);
	  //~ // Wave Speed MINUS
	  //~ double wave_im12 = -u_ij*dsj_x[i][j] - v_ij*dsj_y[i][j] +
	    //~ sound[i][j]*sqrt(dsj_x[i][j]*dsj_x[i][j] + dsj_y[i][j]*dsj_y[i][j]);
	  //~ // Epsilon12 PLUS
	  //~ double eps2pi12 = 0.5*0.25*wave_ip12*
	    //~ std::max(nui[i-1][j], nui[i][j], nui[i+1][j], nui[0][j]);
	  //~ // Epsilon12 MINUS
	  //~ double eps2mi12 = 0.5*0.25*wave_im12*
	    //~ std::max(nui[i-1][j], nui[i][j], nui[i+1][j], nui[i+2][j]);
	//~ }
	//~ if ( i == u[0].size() - 1 ) { // Periodic "i" Fix
	  //~ // Wave speed PLUS
	  //~ double wave_ip12 = u_ij*dsj_x[0][j] + v_ij*dsj_y[0][j] +
	    //~ sound[i][j]*sqrt(dsj_x[0][j]*dsj_x[0][j] + dsj_y[0][j]*dsj_y[0][j]);
	  //~ // Wave Speed MINUS
	  //~ double wave_im12 = -u_ij*dsj_x[i][j] - v_ij*dsj_y[i][j] +
	    //~ sound[i][j]*sqrt(dsj_x[i][j]*dsj_x[i][j] + dsj_y[i][j]*dsj_y[i][j]);
	  //~ // Epsilon12 PLUS
	  //~ double eps2pi12 = 0.5*0.25*wave_ip12*
	    //~ std::max(nui[i-1][j], nui[i][j], nui[0][j], nui[1][j]);
	  //~ // Epsilon12 MINUS
	  //~ double eps2pi12 = 0.5*0.25*wave_im12*
	    //~ std::max(nui[i-1][j], nui[i][j], nui[0][j], nui[1][j]);
	//~ }
	//~ else { // no need for periodic fix
	  //~ // Wave speed PLUS
	  //~ double wave_ip12 = u_ij*dsj_x[i+1][j] + v_ij*dsj_y[i+1][j] +
	    //~ sound[i][j]*sqrt(dsj_x[i+1][j]*dsj_x[i+1][j] + dsj_y[i+1][j]*dsj_y[i+1][j]);
	  //~ // Wave Speed MINUS
	  //~ double wave_im12 = -u_ij*dsj_x[i][j] - v_ij*dsj_y[i][j] +
	    //~ sound[i][j]*sqrt(dsj_x[i][j]*dsj_x[i][j] + dsj_y[i][j]*dsj_y[i][j]);
	  //~ // Epsilon12 PLUS
	  //~ double eps2pi12 = 0.5*0.25*wave_ip12*
	    //~ std::max(nui[i-1][j], nui[i][j], nui[i+1][j], nui[i+2][j]);
	  //~ // Epsilon12 MINUS
	  //~ double eps2mi12 = 0.5*0.25*wave_im12*
	    //~ std::max(nui[i-1][j], nui[i][j], nui[i+1][j], nui[i+2][j]);
	//~ }
	
	//~ // Epsilon4 calc PLUS
	//~ double eps4pi12 = std::max(0, 0.5*0.003906*wave_ip12 - eps2pi12);
	//~ // Epsilon4 calc MINUS
	//~ double eps4mi12 = std::max(0, 0.5*0.003906*wave_im12 - eps2mi12);
	
	//~ // dpi12 calc
	//~ if ( i == 0 ) { // Periodic Fix
	  //~ // DPI12 (plus)
	  //~ for (int k = 0; k < 3; k++) {
	    //~ dpi12[k][i][j] = eps2pi12*(u[k][i+1][j] - u[k][i][j])
	      //~ - eps4pi12*(u[k][i+2][j] - 3.0*u[k][i+1][j]
			  //~ + 3.0*u[k][i][j] - u[k][iSZ][j]);
	  //~ }
	  //~ // DMI12 (minus)
	  //~ for (int k = 0; k < 3; k++) {
	    //~ dmi12[k][i][j] = eps2mi12*(u[k][i+1][j] - u[k][i][j])
	      //~ - eps4mi12*(u[k][i+2][j] - 3.0*u[k][i+1][j]
			  //~ + 3.0*u[k][i][j] - u[k][iSZ][j]);
	  //~ }
	//~ }
	//~ if ( i == u[0].size() - 2 ) { // Periodic Fix
	  //~ // DPI12 (plus)
	  //~ for (int k = 0; k < 3; k++) {
	    //~ dpi12[k][i][j] = eps2pi12*(u[k][i+1][j] - u[k][i][j])
	      //~ - eps4pi12*(u[k][0][j] - 3.0*u[k][i+1][j]
			  //~ + 3.0*u[k][i][j] - u[k][i-1][j]);
	  //~ }
	  //~ // DMI12 (minus)
	  //~ for (int k = 0; k < 3; k++) {
	    //~ dmi12[k][i][j] = eps2mi12*(u[k][i+1][j] - u[k][i][j])
	      //~ - eps4mi12*(u[k][0][j] - 3.0*u[k][i+1][j]
			  //~ + 3.0*u[k][i][j] - u[k][i-1][j]);
	  //~ }
    	//~ }
	//~ if ( i == u[0].size() - 1 ) { // Periodic Fix
	  //~ // DPI12 (plus)
	  //~ for (int k = 0; k < 3; k++) {
	    //~ dpi12[k][i][j] = eps2pi12*(u[k][0][j] - u[k][i][j])
	      //~ - eps4pi12*(u[k][1][j] - 3.0*u[k][0][j]
			  //~ + 3.0*u[k][i][j] - u[k][i-1][j]);
	  //~ }
	  //~ // DMI12 (minus)
	  //~ for (int k = 0; k < 3; k++) {
	    //~ dmi12[k][i][j] = eps2mi12*(u[k][0][j] - u[k][i][j])
	      //~ - eps4mi12*(u[k][1][j] - 3.0*u[k][0][j]
			  //~ + 3.0*u[k][i][j] - u[k][i-1][j]);
	  //~ }
    	//~ }
	//~ else { // no need for periodic fix
	  //~ // DPI12 (plus)
	  //~ for (int k = 0; k < 3; k++) {
	    //~ dpi12[k][i][j] = eps2pi12*(u[k][i+1][j] - u[k][i][j])
	      //~ - eps4pi12*(u[k][i+2][j] - 3.0*u[k][i+1][j]
			  //~ + 3.0*u[k][i][j] - u[k][i-1][j]);
	  //~ }
	  //~ // DMI12 (minus)
	  //~ for (int k = 0; k < 3; k++) {
	    //~ dmi12[k][i][j] = eps2mi12*(u[k][i+1][j] - u[k][i][j])
	      //~ - eps4mi12*(u[k][i+2][j] - 3.0*u[k][i+1][j]
			  //~ + 3.0*u[k][i][j] - u[k][i-1][j]);
	  //~ }
	//~ }
	
	//~ //////////////////////////////////////////////////////////////////////////////
	//~ ///////////////// CALCULATION FOR J'S ////////////////////////////////////////
	//~ ////////////////////////////////////////////// CALCULATION FOR DPJ12, DMJ12 //
	//~ //////////////////////////////////////////////////////////////////////////////
	//~ //////////////////////////////////////////////////////////////////////////////
	//~ /////////////////////////////////// CALCULATION FOR DPJ12 and DMJ12 //////////
	//~ //////////////////////////////////////////////////////////////////////////////
	//~ // Epsilon2 calc
	//~ if ( j == 0 ) { // Double Ghost Fix
	  //~ // Wave Speed PLUS
	  //~ double wave_jp12 = u_ij*dsi_x[i][j+1] + v_ij*dsi_y[i][j+1] +
	    //~ sound[i][j]*sqrt(dsi_x[i][j+1]*dsi_x[i][j+1] + dsi_y[i][j+1]*dsi_y[i][j+1]);
	  //~ // Wave Speed MINUS
	  //~ double wave_jm12 = -u_ij*dsi_x[i][j] - v_ij*dsi_y[i][j] +
	    //~ sound[i][j]*sqrt(dsi_x[i][j]*dsi_x[i][j] + dsi_y[i][j]*dsi_y[i][j]);
	  //~ // Epsilon12 PLUS
	  //~ double eps2pj12 = 0.5*0.25*wave_jp12*
	    //~ std::max(nuj[i][jSZ - 2], nuj[i][j], nuj[i][j+1], nuj[i][j+2]); // CHECK
	  //~ // Epsilon12 MINUS
	  //~ double eps2mj12 = 0.5*0.25*wave_jm12*
	    //~ std::max(nuj[i][jSZ - 2], nuj[i][j], nuj[i][j+1], nuj[i][j+2]); // CHECK
	//~ }
	//~ if ( j == u[0][0].size() - 3 ) { // Double Ghost Fix
	  //~ // Wave speed PLUS
	  //~ double wave_jp12 = u_ij*dsi_x[i][j+1] + v_ij*dsi_y[i][j+1] +
	    //~ sound[i][j]*sqrt(dsi_x[i][j+1]*dsi_x[i][j+1] + dsi_y[i][j+1]*dsi_y[i][j+1]); // CHECK
	  //~ // Wave Speed MINUS
	  //~ double wave_jm12 = -u_ij*dsi_x[i][j] - v_ij*dsi_y[i][j] +
	    //~ sound[i][j]*sqrt(dsi_x[i][j]*dsi_x[i][j] + dsi_y[i][j]*dsi_y[i][j]);
	  //~ // Epsilon12 PLUS
	  //~ double eps2pj12 = 0.5*0.25*wave_jp12*
	    //~ std::max(nuj[i][j - 1], nuj[i][j], nuj[i][j + 1], nuj[i][j + 1]);
	  //~ // Epsilon12 MINUS
	  //~ double eps2mj12 = 0.5*0.25*wave_jm12*
	    //~ std::max(nuj[i][j - 1], nuj[i][j], nuj[i][j + 1], nuj[i][j + 1]);
	//~ }
	//~ if ( j == u[0][0].size() - 4 ) { // Double Ghost Fix
	  //~ // Wave speed PLUS
	  //~ double wave_jp12 = u_ij*dsi_x[i][j+1] + v_ij*dsi_y[i][j+1] +
	    //~ sound[i][j]*sqrt(dsi_x[i][j+1]*dsi_x[i][j+1] + dsi_y[i][j+1]*dsi_y[i][j+1]);
	  //~ // Wave Speed MINUS
	  //~ double wave_jm12 = -u_ij*dsi_x[i][j] - v_ij*dsi_y[i][j] +
	    //~ sound[i][j]*sqrt(dsi_x[i][j]*dsi_x[i][j] + dsi_y[i][j]*dsi_y[i][j]);
	  //~ // Epsilon12 PLUS
	  //~ double eps2pj12 = 0.5*0.25*wave_jp12*
	    //~ std::max(nuj[i][j - 1], nuj[i][j], nuj[i][j + 1], nuj[1][j + 2]);
	  //~ // Epsilon12 MINUS
	  //~ double eps2pj12 = 0.5*0.25*wave_jm12*
	    //~ std::max(nuj[i][j - 1], nuj[i][j], nuj[i][j + 1], nuj[1][j + 2]);
	//~ }
	//~ else { // no need for ghost fix
	  //~ // Wave speed PLUS
	  //~ double wave_jp12 = u_ij*dsi_x[i][j+1] + v_ij*dsi_y[i][j+1] +
	    //~ sound[i][j]*sqrt(dsi_x[i][j+1]*dsi_x[i][j+1] + dsi_y[i][j+1]*dsi_y[i][j+1]);
	  //~ // Wave Speed MINUS
	  //~ double wave_jm12 = -u_ij*dsj_x[i][j] - v_ij*dsj_y[i][j] +
	    //~ sound[i][j]*sqrt(dsj_x[i][j]*dsj_x[i][j] + dsj_y[i][j]dsj_y[i][j]);
	  //~ // Epsilon12 PLUS
	  //~ double eps2pj12 = 0.5*0.25*wave_jp12*
	    //~ std::max(nuj[i][j - 1], nuj[i][j], nuj[i][j+1], nuj[i][j+2]);
	  //~ // Epsilon12 MINUS
	  //~ double eps2mj12 = 0.5*0.25*wave_jm12*
	    //~ std::max(nuj[i][j - 1], nuj[i][j], nuj[i][j+1], nuj[i][j+2]);
	//~ }
	
	//~ // Epsilon4 calc PLUS
	//~ double eps4pj12 = std::max(0, 0.5*0.003906*wave_jp12 - eps2pj12);
	//~ // Epsilon4 calc MINUS
	//~ double eps4mj12 = std::max(0, 0.5*0.003906*wave_jm12 - eps2mj12);
	
	//~ // dpj12 calc
	//~ if ( j == 0 ) { // Ghost Fix
	  //~ // DPJ12 (plus)
	  //~ for (int k = 0; k < 3; k++) {
	    //~ dpj12[k][i][j] = eps2pj12*(u[k][i][j+1] - u[k][i][j])
	      //~ - eps4pj12*(u[k][i][j+2] - 3.0*u[k][i][j+1]
			  //~ + 3.0*u[k][i][j] - u[k][i][j-1]);
	  //~ }
	  //~ // DMJ12 (minus)
	  //~ for (int k = 0; k < 3; k++) {
	    //~ dmj12[k][i][j] = eps2mj12*(u[k][i][j+1] - u[k][i][j])
	      //~ - eps4mj12*(u[k][i][j+2] - 3.0*u[k][i][j+1]
			  //~ + 3.0*u[k][i][j] - u[k][i][j-1]);
	  //~ }
	//~ }
	//~ if ( j == u[0][0].size() - 4 ) { // Ghost Fix
	  //~ // DPJ12 (plus)
	  //~ for (int k = 0; k < 3; k++) {
	    //~ dpj12[k][i][j] = eps2pj12*(u[k][i][j+1] - u[k][i][j])
	      //~ - eps4pj12*(u[k][i][j+2] - 3.0*u[k][i][j+1]
			  //~ + 3.0*u[k][i][j] - u[k][i][j-1]); // CHECK
	  //~ }
	  //~ // DMJ12 (minus)
	  //~ for (int k = 0; k < 3; k++) {
	    //~ dmj12[k][i][j] = eps2mj12*(u[k][i][j+1] - u[k][i][j])
	      //~ - eps4mj12*(u[k][i][j+2] - 3.0*u[k][i][j+1]
			  //~ + 3.0*u[k][i][j] - u[k][i][j-1]); // CHECK
	  //~ }
    	//~ }
	//~ if ( j == u[0][0].size() - 3 ) { // Ghost Fix
	  //~ // DPJ12 (plus)
	  //~ for (int k = 0; k < 3; k++) {
	    //~ dpj12[k][i][j] = eps2pj12*(u[k][i][j+1] - u[k][i][j])
	      //~ - eps4pj12*(u[k][i][j+2] - 3.0*u[k][i][j+1]
			  //~ + 3.0*u[k][i][j] - u[k][i][j-1]); // CHECK
	  //~ }
	  //~ // DMJ12 (minus)
	  //~ for (int k = 0; k < 3; k++) {
	    //~ dmj12[k][i][j] = eps2mj12*(u[k][i][j+1] - u[k][i][j])
	      //~ - eps4mj12*(u[k][i][j+2] - 3.0*u[k][i][j+1]
			  //~ + 3.0*u[k][i][j] - u[k][i][j-1]); // CHECK
	  //~ }
    	//~ }
	//~ else { // no need for ghost fix
	  //~ // DPJ12 (plus)
	  //~ for (int k = 0; k < 3; k++) {
	    //~ dpj12[k][i][j] = eps2pj12*(u[k][i][j+1] - u[k][i][j])
	      //~ - eps4pj12*(u[k][i][j+2] - 3.0*u[k][i][j+1]
			  //~ + 3.0*u[k][i][j] - u[k][i][j-1]);
	  //~ }
	  //~ // DMJ12 (minus)
	  //~ for (int k = 0; k < 3; k++) {
	    //~ dmj12[k][i][j] = eps2mj12*(u[k][i][j+1] - u[k][i][j])
	      //~ - eps4mj12*(u[k][i][j+2] - 3.0*u[k][i][j+1]
			  //~ + 3.0*u[k][i][j] - u[k][i][j-1]);
	  //~ }
	//~ }
	
    //~ } // Close Loop over J
    //~ } // Close Loop over I
    
    
    // Calculate residuals
    // First calculate wave-speeds and normalized pressure gradients (???)
    
    for (int i=0; i < u[0].size(); i++){
        int im = i - 1;
        int ip = i + 1;
        if (i == 0) {
            // Periodic boundary (i-1 = i_max)
            im = u[0].size() - 1;
        } else if (i == u[0].size() - 1) {
            // Periodic boundary (i-1 = i_max)
            ip = 0;
        }
        
        for (int j=0; j < u[0][0].size()-2; j++) {
            int jm = j - 1;
            int jp = j + 1;
            if (j == 0) {
                // Call airfoil ghost cell 
                jm = u[0][0].size() - 2;
            } else if (j == u[0][0].size() - 3) {
                // Call external boundary ghost cell
                jp = u[0][0].size() - 1;
            }
            
            for (int var=0; var<u.size(); var++) {
                // average fluxes at each wall between cells
                // then dot with the wall normal
                double Fmi12 = 0.5*(F[var][i][j] + F[var][im][j])*dsj_x[i][j];
                double Fpi12 = 0.5*(F[var][i][j] + F[var][ip][j])*dsj_x[i+1][j];
                double Fmj12 = 0.5*(F[var][i][j] + F[var][i][jm])*dsi_x[i][j];
                double Fpj12 = 0.5*(F[var][i][j] + F[var][i][jp])*dsi_x[i][j+1];
                                                                 
                double Gmi12 = 0.5*(G[var][i][j] + G[var][im][j])*dsj_y[i][j];
                double Gpi12 = 0.5*(G[var][i][j] + G[var][ip][j])*dsj_y[i+1][j];
                double Gmj12 = 0.5*(G[var][i][j] + G[var][i][jm])*dsi_y[i][j];
                double Gpj12 = 0.5*(G[var][i][j] + G[var][i][jp])*dsi_y[i][j+1];
                
                //~ if (j==u[0][0].size()-3 && var==1){
                    //~ std::cout << "i = " << i << std::endl;
                    //~ std::cout << "Fmi12 = " << Fmi12 << std::endl;
                    //~ std::cout << "Fpi12 = " << Fpi12 << std::endl;
                    //~ std::cout << "Fmj12 = " << Fmj12 << std::endl;
                    //~ std::cout << "Fpj12 = " << Fpj12 << std::endl;
                    //~ std::cout << "Gmi12 = " << Gmi12 << std::endl;
                    //~ std::cout << "Gpi12 = " << Gpi12 << std::endl;
                    //~ std::cout << "Gmj12 = " << Gmj12 << std::endl;
                    //~ std::cout << "Gpj12  = " << Gpj12 << std::endl;
                //~ }
                // Air foil boundary flux condition
                if (j==0) {
                    if (var == 0 || var == 3) {
                        // the flux of rho and E normal to air foil = 0
                        r[var][i][j]  = -1.*(Fmi12 + Gmi12) ; //+ dmi12[var][i][j]; 
                        r[var][i][j] +=      Fpi12 + Gpi12  ; //- dpi12[var][i][j];
                        r[var][i][j] += 0.;                 ; //
                        r[var][i][j] +=      Fpj12 + Gpj12  ; //- dpj12[var][i][j];
                    } else if (var == 1) {                   //
                        // the flux of momentum = p         ; //
                        r[var][i][j]  = -1.*(Fmi12 + Gmi12) ; //+ dmi12[var][i][j]; 
                        r[var][i][j] +=      Fpi12 + Gpi12  ; //- dpi12[var][i][j];
                        r[var][i][j] += -1.*p[i][j]*dsi_x[i][j]; //
                        r[var][i][j] +=      Fpj12 + Gpj12  ; //- dpj12[var][i][j];
                    } else if (var == 2){                    //
                        // the flux of momentum = p         ; //
                        r[var][i][j]  = -1.*(Fmi12 + Gmi12) ; //+ dmi12[var][i][j]; 
                        r[var][i][j] +=      Fpi12 + Gpi12  ; //- dpi12[var][i][j];
                        r[var][i][j] += -1.*p[i][j]*dsi_y[i][j];; //
                        r[var][i][j] +=      Fpj12 + Gpj12  ; //- dpj12[var][i][j];
                    }
                } else {
                    // Sum up fluxes normally 
                    r[var][i][j]  = -1.*(Fmi12 + Gmi12) ; //+ dmi12[var][i][j]; 
                    r[var][i][j] +=      Fpi12 + Gpi12  ; //- dpi12[var][i][j];
                    r[var][i][j] += -1.*(Fmj12 + Gmj12) ; //+ dmj12[var][i][j];
                    r[var][i][j] +=      Fpj12 + Gpj12  ; //- dpj12[var][i][j];
                }
            }
        }
    }
    //~ writeOutput(r);
}

void calcTau(const std::vector<std::vector<std::vector<double> > > &u, 
	     const int i, const int j, double &tau_ij,
	     std::vector<std::vector<double> > &dsi_x,
	     std::vector<std::vector<double> > &dsi_y,
         std::vector<std::vector<double> > &dsj_x,
         std::vector<std::vector<double> > &dsj_y)
{
    // Calculate local time step
    // Given parameters
    double cfl = 0.5;
    double gamma = 1.4;
    
    // Find local u, v and c
    double rho_ij = u[0][i][j];
    double u_ij = u[1][i][j]/rho_ij;
    double v_ij = u[2][i][j]/rho_ij;
    double p_ij = (gamma - 1)*(u[3][i][j] - 0.5*rho_ij*(u_ij*u_ij + v_ij*v_ij));
    double c_ij = std::sqrt(gamma*p_ij/rho_ij);
    
    //~ std::cout << "P_ij = " << p_ij << std::endl;
    // Average wall normals
    double dsi_x_avg = 0.5*(dsi_x[i][j] + dsi_x[i][j+1]);
    double dsi_y_avg = 0.5*(dsi_y[i][j] + dsi_y[i][j+1]);
    double dsj_x_avg = 0.5*(dsj_x[i][j] + dsj_x[i+1][j]);
    double dsj_y_avg = 0.5*(dsj_y[i][j] + dsj_y[i+1][j]);
    
    // Calculate max possible speed
    double uc = u_ij + c_ij;
    double vc = v_ij + c_ij;
    
    // Calculate time step
    tau_ij = cfl/(std::abs(dsi_x_avg*uc + dsi_y_avg*vc) 
                                + std::abs(dsj_x_avg*uc + dsj_y_avg*vc));
    
}

void tempInt(const std::vector<std::vector<std::vector<double> > > &u,
            const std::vector<std::vector<std::vector<double> > > &r,
            std::vector<std::vector<double> > &dsi_x,
            std::vector<std::vector<double> > &dsi_y,
            std::vector<std::vector<double> > &dsj_x,
            std::vector<std::vector<double> > &dsj_y,
            const double alpha,
	    std::vector<std::vector<std::vector<double> > > &u_new)
{
    // Computes RK step given variable vector u, residual from spacial 
    // integrator r, a local time step (calculated in calcTau) and RK step weight alpha.
    
    double tau_ij = 0.;
    for (int i=0; i<u[0].size(); i++){
        for (int j=0; j<u[0][0].size()-2; j++) {
            // Calculate local time step
	    calcTau(u, i, j, tau_ij, dsi_x, dsi_y, dsj_x, dsj_y);
            
            // Calc RK step
            u_new[0][i][j] = u[0][i][j] - alpha*tau_ij*r[0][i][j];
            u_new[1][i][j] = u[1][i][j] - alpha*tau_ij*r[1][i][j];
            u_new[2][i][j] = u[2][i][j] - alpha*tau_ij*r[2][i][j];
            u_new[3][i][j] = u[3][i][j] - alpha*tau_ij*r[3][i][j];
            //~ std::cout << tau_ij << std::endl;
        }
    }
}

// Function to set exterior BCs using Riemann invariants
void setExteriorBC(std::vector< std::vector< std::vector<double> > > &u,
		   std::vector<std::vector<double> > &dsi_x,
		   std::vector<std::vector<double> > &dsi_y,
		   std::vector<std::vector<double> > &dsj_x,
		   std::vector<std::vector<double> > &dsj_y,
		   const int &N_col,
		   const double &u_0,
		   const double &c_0,
		   const double &rho_0,
		   const double &p_0,
		   const double &gamma) {
  
  ///// Loop over the upper row of grid
  // Initialize stuff
  double uin = 0.; // u dot n interior
  double u0n = 0.; // u dot n infty
  //double uit = 0.; // u tangential interior
  double u0t = 0.; // u tangential infty
  double uvel = 0.; // ux interior
  double vvel = 0.; // uy interior
  double norm = 0.; // For normalizing normal vector
  double nx = 0.;
  double ny = 0.;
  double R_inf = 0.; // Riemann infty
  double R_int = 0.; // Riemann interior
  double ci = 0.; // interior sound speed
  double cb = 0.; // boundary sound speed
  double unb = 0.; // u dot n boundary
  double rhob = 0.; // density boundary
  double pb = 0.; // pressure boundary
  double pi = 0.; // pressure interior
  int I_max = u[0].size();
  int J_max = u[0][0].size();
  
  for (int i = 0; i < I_max; i++) {
    // Normalize "j" unit vectors
    norm = sqrt(dsi_x[i][J_max - 3]*dsi_x[i][J_max - 3]
		+ dsi_y[i][J_max - 3]*dsi_y[i][J_max - 3]);
    
    nx = -1.0*dsi_x[i][J_max - 3]/norm; // Reverse for convention
    ny = -1.0*dsi_y[i][J_max - 3]/norm; // Reverse for convention
    
    // Check for inflow or outflow
    uvel = u[1][i][J_max - 3]/u[0][i][J_max - 3];
    vvel = u[2][i][J_max - 3]/u[0][i][J_max - 3];
    
    // Compute u dot n (normal component of velocity)
    uin = uvel*nx + vvel*ny; // interior
    u0n = u_0*nx; // infty
    
    // Compute u tangential
    //uit = uvel*ny - vvel*nx; // interior
    u0t = u_0*ny; // infty
    
    // Compute interior speed of sound
    ci = sqrt(gamma*(gamma - 1)*(
	      u[3][i][J_max - 3]/u[0][i][J_max - 3]
	      - 0.5*(uvel*uvel + vvel*vvel)));
    
    // Compute interior pressure
    pi = u[0][i][J_max - 3]*ci*ci/gamma;

    if (i == 120) {
      std::cout << "uin = " << uin << std::endl;
    }
    
    // Continue as inflow or outflow
    if (uin > 0) { // inflow BC
      // Compute R_inf (infty, plus)
      // and R_int (interior, minus)
      R_inf = abs(u0n) + 2.0/(gamma - 1.0)*c_0;
      R_int = abs(uin) - 2.0/(gamma - 1.0)*ci;
      
      // Compute Primitive Boundary Values
      unb = 0.5*(R_inf + R_int);
      cb = 0.25*(gamma - 1.0)*(R_inf - R_int);
      rhob = pow(cb*cb*pow(rho_0, gamma)/(gamma*p_0), 1./(gamma - 1.0));
      pb = rhob*cb*cb/gamma;
      
      // Inflow
      if (i < 64 && i > 32) {
	
	// Compute Conserved Boundary Values
	// density
	u[0][i][J_max - 1] = rhob;
	// x velocity
	u[1][i][J_max - 1] = rhob*(unb*nx + u0t*ny);
	// y velocity
	u[2][i][J_max - 1] = rhob*(unb*ny - u0t*nx);
	// Energy
	u[3][i][J_max - 1] = rhob*((1./(gamma - 1.0))*(pb/rhob)
			        + 0.5*(unb*unb + u0t*u0t));
      }
      if (i >= 64) {
	
	// Compute Conserved Boundary Values
	// density
	u[0][i][J_max - 1] = rhob;
	// x velocity
	u[1][i][J_max - 1] = rhob*(unb*nx + u0t*ny);
	// y velocity
	u[2][i][J_max - 1] = rhob*(unb*ny - u0t*nx);
	// Energy
	u[3][i][J_max - 1] = rhob*((1./(gamma - 1.0))*(pb/rhob)
				   + 0.5*(unb*unb + u0t*u0t));
      } 
    }
    
    if (uin < 0) { // outflow BC
      // Compute R_inf (infty, minus)
      // and R_int (interior, plus)
      R_inf = -1.0*abs(u0n) - 2.0/(gamma - 1.0)*c_0;
      R_int = -1.0*abs(uin) + 2.0/(gamma - 1.0)*ci;
      
      // Compute Primitive Boundary Values
      unb = 0.5*(R_inf + R_int);
      cb = 0.25*(gamma - 1.0)*(R_inf - R_int);
      rhob = pow(cb*cb*pow(rho_0, gamma)/(gamma*p_0), 1./(gamma - 1.0));
      pb = rhob*cb*cb/gamma;
      
      // Outflow
      if (i <= 32) {
	
	// Compute Conserved Boundary Values
	// density
	u[0][i][J_max - 1] = rhob;
	// x velocity
	u[1][i][J_max - 1] = rhob*(unb*nx + u0t*ny);
	// y velocity
	u[2][i][J_max - 1] = rhob*(unb*ny - u0t*nx);
	//u[2][i][J_max - 1] = 0;
	// Energy
	u[3][i][J_max - 1] = rhob*((1./(gamma - 1.0))*(pb/rhob)
				   + 0.5*(unb*unb + u0t*u0t));
      }
      if (i > 96) {
	
	// Compute Conserved Boundary Values
	// density
	u[0][i][J_max - 1] = rhob;
	// x velocity
	u[1][i][J_max - 1] = rhob*(unb*nx + u0t*ny);
	// y velocity
	u[2][i][J_max - 1] = rhob*(unb*ny - u0t*nx);
	//u[2][i][J_max - 1] = 0;
	// Energy
	u[3][i][J_max - 1] = rhob*((1./(gamma - 1.0))*(pb/rhob)
				   + 0.5*(unb*unb + u0t*u0t));
      }
    }
  }
}

// Perscribe values in airfoil ghost cells
void setAirfoilBC ( std::vector<std::vector<std::vector<double> > > &u,
                    const std::vector<std::vector<double> > &dsi_x,
                    const std::vector<std::vector<double> > &dsi_y)
{
    // Copy interior cell values, but setting normal velocity to zero.
    int b = u[0][0].size() - 2;  // boundary cell

    for (int i=0; i<u[0].size(); i++) {
        // Calculate unit normals
        double s_mag = std::sqrt(dsi_x[i][0]*dsi_x[i][0] + dsi_y[i][0]*dsi_y[i][0]);
        double nx = dsi_x[i][0] / s_mag;
        double ny = dsi_y[i][0] / s_mag;
        
        // Find normal velocity
        double u_norm = u[1][i][b-1] / u[0][i][b-1] * nx;
        double v_norm = u[2][i][b-1] / u[0][i][b-1] * ny;
        
        // Find tangential velocity
        double u_tan =  u[1][i][b-1] / u[0][i][b-1] * ny;
        double v_tan = -u[2][i][b-1] / u[0][i][b-1] * nx;
        
        // Swap sign of normal
        u_norm = -u_norm;
        v_norm = -v_norm;
        
        // Assign ghost cell values
        u[0][i][b] = u[0][i][b-1];  // same as interior
        u[1][i][b] = u[0][i][b-1]*(u_norm + u_tan);  // same tan, opps normal
        u[2][i][b] = u[0][i][b-1]*(v_norm + v_tan);  // same tan, ops normal
        u[3][i][b] = u[3][i][b-1];  // same as interior
    }
}


int main()
{
    // Run parameters
    double cfl      = 0.5;
    double R_min    = 1.0e-6; // Target residual value
    double vars     = 4; // Number of conserved variables
    
    // Free stream values (used in ICs and BCs)
    double gamma    = 1.4;
    double rho_0    = 0.414;
    double u_0      = 255.;
    double v_0      = 0.;
    double c_0      = 300.;
    //double E_0      = 193.;
    double p_0      = 26.5e3;
    
    std::cout << "p_0" << p_0 << std::endl;
    
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
    std::vector< std::vector<double> > u_temp;
    std::vector< std::vector< std::vector<double> > > u;
    
    ////// Size the vectors
    std::vector<double> col (N_col); // For vectors that hold verticie values
    std::vector<double> colm (N_col-1); // For vectors that hold cell centered values
    std::vector<double> colp (N_col+1); // For vectors with ghost cells (just u)
    grid_x.push_back(col);
    grid_y.push_back(col);
    dsj_x.push_back(colm);
    dsj_y.push_back(colm);
    
    for (int row = 0; row < N_row - 1; row++) {
        grid_x.push_back(col);
        grid_y.push_back(col);
        x.push_back(colm);
        y.push_back(colm);
        omega.push_back(colm);
        dsi_x.push_back(col);
        dsi_y.push_back(col);
        dsj_x.push_back(colm);
        dsj_y.push_back(colm);
        u_temp.push_back(colp);
    }
    
    // create u vector (N_col+1 -> airfoil ghost cell, N_col+2 -> exterior ghost cell)
    for (int var = 0; var<vars; var++){
        u.push_back(u_temp);
    }
    
    std::cout << "u[].size() = " << u[0].size() << " u[][].size = " << u[0][0].size() << std::endl; 
    std::cout << "grid_x.size() = " << grid_x.size() << " grid_x[].size = " << grid_x[0].size() << std::endl;
    std::cout << "x.size() = " << x.size() << " x[].size = " << x[0].size() << std::endl;
    std::cout << "dsi_x.size() = " << dsi_x.size() << " dsi_x[].size = " << dsi_x[0].size() << std::endl; 
    std::cout << "dsj_x.size() = " << dsj_x.size() << " dsj_x[].size = " << dsj_x[0].size() << std::endl; 

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
    writeCellNormals(dsi_x, dsi_y, dsj_x, dsj_y);
    
    // Set up initial conditions
    setIC(u, rho_0, u_0, v_0, p_0, gamma);
    //~ writeOutput(u);
    
    // Run sim (Using standard RK4 in time and Jameson artifical viscosty in space)
    double R_max = 1.0;
    double R_avg = 1.0;
    
    //~ while (R_max > 1.0e-6) {
    for (int i=0; i<500; i++) {
        // Create copies of u for RK steps, plus an r to store residuals;
        std::vector< std::vector< std::vector<double> > > u_1(u);
        std::vector< std::vector< std::vector<double> > > u_2(u);
        std::vector< std::vector< std::vector<double> > > u_3(u);
        std::vector< std::vector< std::vector<double> > > u_4(u);
        std::vector< std::vector< std::vector<double> > > r(u);
        
        setExteriorBC(u, dsi_x, dsi_y, dsj_x, dsj_y, N_col, u_0, c_0, rho_0, p_0, gamma);
        setAirfoilBC(u, dsi_x, dsi_y);
	
	writeOutput(u);
        
        // First RK step
        // Calculate residual
        spaceInt(u, omega, dsi_x, dsi_y, dsj_x, dsj_y, r, gamma);
        //~ writeOutput(r);
        // Calculate first RK step
        double alpha = 0.25;
        tempInt(u, r, dsi_x, dsi_y, dsj_x, dsj_y, alpha, u_1);
        
        setAirfoilBC(u_1, dsi_x, dsi_y);
        setExteriorBC(u_1, dsi_x, dsi_y, dsj_x, dsj_y, N_col, u_0, c_0, rho_0, p_0, gamma);
        
        writeOutput(u_1);
        
        // Second RK step
        // Calculate residual using u_1
        //spaceInt(u_1, omega, dsi_x, dsi_y, dsj_x, dsj_y, r, gamma); 
        // Calculate RK step
        alpha = 1./3.;
        
        tempInt(u, r, dsi_x, dsi_y, dsj_x, dsj_y, alpha, u_2);
        
        setExteriorBC(u_2, dsi_x, dsi_y, dsj_x, dsj_y, N_col, u_0, c_0, rho_0, p_0, gamma);
        setAirfoilBC(u_2, dsi_x, dsi_y);
	
        // Third RK step
        // Calculate residual using u_2
        spaceInt(u_2, omega, dsi_x, dsi_y, dsj_x, dsj_y, r, gamma);
        // Calculate RK step
        alpha = 0.5;
        tempInt(u, r, dsi_x, dsi_y, dsj_x, dsj_y, alpha, u_3);
        
        setExteriorBC(u_3, dsi_x, dsi_y, dsj_x, dsj_y, N_col, u_0, c_0, rho_0, p_0, gamma);
        setAirfoilBC(u_3, dsi_x, dsi_y);
        
        // Fourth RK step
        // Calculate residual using u_3
        spaceInt(u_3, omega, dsi_x, dsi_y, dsj_x, dsj_y, r, gamma);
        // Calculate RK step
        alpha = 1.;
        tempInt(u, r, dsi_x, dsi_y, dsj_x, dsj_y, alpha, u_4);
        
        setExteriorBC(u_4, dsi_x, dsi_y, dsj_x, dsj_y, N_col, u_0, c_0, rho_0, p_0, gamma);
        setAirfoilBC(u_4, dsi_x, dsi_y);
        
        
        // Update original vector u
        u = u_4;
        
        // Calculate total residual
        R_max = 0.0;
        R_avg = 0.0;
        for (int i=0; i<r[0].size(); i++) {
            for (int j=0; j<r[0][0].size()-2; j++) {
                double r_max = 0.0;
                if (r[0][i][j] > r[1][i][j] || r[2][i][j] > r[3][i][j]) {
                    r_max = std::max(r[0][i][j], r[2][i][j]);
                } else {
                    r_max = std::max(r[1][i][j], r[3][i][j]);
                }
                if (r_max > R_max) {
                    R_max = r_max;
                }
                R_avg = R_avg + R_max;
            }
        }
        R_avg = R_avg / (r[0].size() * r[0][0].size());
        std::cout << "R_max = " << R_max << " R_avg = " << R_avg << std::endl;
    }
    
    writeOutput(u);
    
    return 0;
}
