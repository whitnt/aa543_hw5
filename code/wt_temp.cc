            F[0][i][j] = u_1;
            F[1][i][j] = u_1*u_1/u_0 + p;
            F[2][i][j] = u_1*u_2/u_0;
            F[3][i][j] = u_1*(u_3 + p)/u_0;
            G[0][i][j] = u_2;
            G[1][i][j] = u_1*u_2/u_0;
            G[2][i][j] = u_2*u_2/u_0 + p;
            G[3][i][j] = u_2*(u_3 + p)/u_0;




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
            
            for (int var; var<u.size(); var++) {

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
                
                // Air foil boundary flux condition
                if (i==0) {
                    if (var == 0 || var == 3) {
                        // the flux of rho and E normal to air foil = 0
                        r[var][i][j]  = -1.*(Fmi12 + Gmi12) - dmi12[var][i][j]; 
                        r[var][i][j] +=      Fpi12 + Gpi12  - dpi12[var][i][j];
                        r[var][i][j] += 0.;
                        r[var][i][j] +=      Fpj12 + Gpj12  - dpj12[var][i][j];
                    } else if (var == 1) {
                        // the flux of momentum = p
                        r[var][i][j]  = -1.*(Fmi12 + Gmi12) - dmi12[var][i][j]; 
                        r[var][i][j] +=      Fpi12 + Gpi12  - dpi12[var][i][j];
                        r[var][i][j] += p[i][j]*dsi_x[i][j];
                        r[var][i][j] +=      Fpj12 + Gpj12  - dpj12[var][i][j];
                    } else if (var == 2){
                        // the flux of momentum = p
                        r[var][i][j]  = -1.*(Fmi12 + Gmi12) - dmi12[var][i][j]; 
                        r[var][i][j] +=      Fpi12 + Gpi12  - dpi12[var][i][j];
                        r[var][i][j] += p[i][j]*dsi_y[i][j];
                        r[var][i][j] +=      Fpj12 + Gpj12  - dpj12[var][i][j];
                    }
                } else if {
                    // Sum up fluxes normally 
                    r[var][i][j]  = -1.*(Fmi12 + Gmi12) - dmi12[var][i][j]; 
                    r[var][i][j] +=      Fpi12 + Gpi12  - dpi12[var][i][j];
                    r[var][i][j] += -1.*(Fmj12 + Gmj12) - dmj12[var][i][j];
                    r[var][i][j] +=      Fpj12 + Gpj12  - dpj12[var][i][j];
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
        double s_mag = std::sqrt(dsi_x[i][0]*dsi_x[i][0] + dsi_y[i][0]*dsi_y[i][0]) 
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

// Calculate total residual
void calcTotalResidual (const std::vector<std::vector<std::vector<double> > > &u)
{
    // Simple average
    
    // Weighted average
    
    // Max value
    for (int i=0; i<r[0].size(); i++) {
        for (int j=0;
    }
}
