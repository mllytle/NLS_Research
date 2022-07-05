function y = NLSE1d_msd(phi,params)

    %% Get params 
    m = params.m;
%     Vfunct = params.V;
%     V = Vfunct(phi);
    V = params.Vstat;
    mu = params.mu;
    h = params.h;
%     a = params.a;
%     s = params.s;

    % Pre-allocate nonlinear residual and solution:
    y_r = zeros(m,1);
    y_i = zeros(m,1);

    %% real and imag field parts
    phi_r = phi(1:m);
    phi_i = phi(m+1:2*m);

    %% Density and common term
    dens = phi_r.^2 + phi_i.^2;
    comm =  dens + V' - mu*ones(m,1);
    
    %% Compute the 1D Laplacians inside the domain 
    % (i.e., from j=1 to the Nth grid point):
    d2Xdx2 = diff(phi_r,2) / h^2; 
    d2Ydx2 = diff(phi_i,2) / h^2;

    %% Compute the common term 
    % (see, the term in the square brackets in Eq. (3.4) in Hermano Ricardo's paper)--I am calling it \Omega:
    term_l = ( d2Xdx2(1).*phi_r(2)+d2Ydx2(1).*phi_i(2) ) / dens(2);
    term_r = ( d2Xdx2(m-2).*phi_r(m-1)+d2Ydx2(m-2).*phi_i(m-1) ) / dens(m-1);
    Omega_l = term_l - 2 * ( dens(2) - dens(1) + V(2) - V(1) );
    Omega_r = term_r - 2 * ( dens(m-1) - dens(m) + V(m-1) - V(m) );

    %% Compute the 2nd-order derivatives of real, imag, and endpoints
    Xdd_l = Omega_l * phi_r(1); 
    Xdd_r = Omega_r * phi_r(m);
    Ydd_l = Omega_l * phi_i(1); 
    Ydd_r = Omega_r * phi_i(m);

    d2Xdx2 = [ Xdd_l; d2Xdx2; Xdd_r ];
    d2Ydx2 = [ Ydd_l; d2Ydx2; Ydd_r ];

    %% Return the system of nonlinear equations
    y_r = -0.5 * d2Xdx2 + comm.*phi_r;
    y_i = -0.5 * d2Ydx2 + comm.*phi_i;

    y = [y_r; y_i];

end