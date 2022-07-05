function jac = jac_32(phi, params)
% x is concatenation of R and Imag comps of wave eq. at disc. points
    
    %% get params
%     s = params.s;
    m = params.m;
%     a = params.a;
%     Vfunct = params.V;
%     V = Vfunct(params.xdom); %not phi?
    V = params.Vstat; %***
    mu = params.mu;
    h = params.h;

    phi_r = phi(1:m); %X
    phi_rr2 = phi_r.^2; %X2
    phi_i = phi(m+1:2*m); %Y
    phi_ii2 = phi_i.^2; %Y2

    A11 = zeros(m);
    A12 = zeros(m);
    A21 = zeros(m);
    A22 = zeros(m);

    %% Get vars
    phi_r0 = phi_r(1); %XX0
    phi_r1 = phi_r(2);
    phi_r2 = phi_r(3);
    phi_rnpm = phi_r(m-2);
    phi_rn = phi_r(m-1);
    phi_rnpp = phi_r(m); %XXNp1

    phi_i0 = phi_i(1); %YY0
    phi_i1 = phi_i(2);
    phi_i2 = phi_i(3);
    phi_inpm = phi_i(m-2);
    phi_in = phi_i(m-1);
    phi_inpp = phi_i(m); %YYNp1

    h2 = h^2; %dx2

    V0 = V(1);
    V1 = V(2);
    Vn = V(m-1); %VN
    Vnpp = V(m); %Vnpp

    comm =  V - mu*ones(m,1); %***

    % Derivatives: (to construct leading and trailing rows)

    %% I eq w/t respect to R
    Omega_l = ((phi_r2 - 2*phi_r1 + phi_r0)*phi_r1 + (phi_i2 - 2*phi_i1 + phi_i0)*phi_i1)...
        /(h2*(phi_r1^2 + phi_i1^2)) - 2*(phi_r1^2 + phi_i1^2 - phi_r0^2 - phi_i0^2 + V1 - V0);
    dOmega_l_X0 = 4*phi_r0 + phi_r1/(h2*(phi_r1^2 + phi_i1^2));
    dOmega_l_X1 = -4*phi_r1 + (phi_r0 - 4*phi_r1 + phi_r2)/(h2*(phi_r1^2 + phi_i1^2)) ...
        - ( 2*phi_r1*(phi_r1*(phi_r0 - 2*phi_r1 + phi_r2) + phi_i1*(phi_i0 - 2*phi_i1 + phi_i2)))/(h2*(phi_r1^2 + phi_i1^2)^2);
    dOmega_l_X2 = phi_r1/(h2*(phi_r1^2 + phi_i1^2));

    %% I eq w/t respect to I
    dOmega_l_Y0 = 4*phi_i0 + phi_i1/(h2*(phi_r1^2 + phi_i1^2));
    dOmega_l_Y1 = -4*phi_i1 + (phi_i0 - 4*phi_i1 + phi_i2)/(h2*(phi_r1^2 + phi_i1^2)) ...
        - (2*phi_i1*(phi_r1*(phi_r0 - 2*phi_r1 + phi_r2) + phi_i1*(phi_i0 - 2*phi_i1 + phi_i2)))...
        /(h2*(phi_r1^2 + phi_i1^2)^2); 
     dOmega_l_Y2 = phi_i1/(h2*(phi_r1^2 + phi_i1^2)); 
    
    %% R eq w/t respect to R
    Omega_r = ((phi_rnpp - 2*phi_rn + phi_rnpm)*phi_rn + (phi_inpp - 2*phi_in + phi_inpm)*phi_in) ...
        /(h2*(phi_rn^2 + phi_in^2)) - 2*(phi_rn^2 + phi_in^2 - phi_rnpp^2 - phi_inpp^2 + Vn - Vnpp);
    dOmega_r_XNm1 = phi_rn/(h2*(phi_rn^2 + phi_in^2)); 
    dOmega_r_XN = -4*phi_rn + (-4*phi_rn + phi_rnpm + phi_rnpp)/(h2*(phi_rn^2 + phi_in^2))...
        - (2*phi_rn*(phi_rn*(-2*phi_rn + phi_rnpm + phi_rnpp) ...
        + phi_in * (-2*phi_in + phi_inpm + phi_inpp)))...
        /(h2*(phi_rn^2 + phi_in^2)^2);
    dOmega_r_XNp1 = 4*phi_rnpp + phi_rn/(h2*(phi_rn^2 + phi_in^2));

    %% R eq w/t respect to I
    dOmega_r_YNm1 = phi_in/(h2*(phi_rn^2 + phi_in^2)); 
    dOmega_r_YN = -4*phi_in + (-4*phi_in + phi_inpm + phi_inpp)/(h2*(phi_rn^2 + phi_in^2)) ...
        - (2*phi_in*(phi_rn*(-2*phi_rn + phi_rnpm + phi_rnpp) ...
        + phi_in*(-2*phi_in + phi_inpm + phi_inpp)))/(h2*(phi_rn^2 + phi_in^2)^2);
    dOmega_r_YNp1 = phi_in/(h2*(phi_rn^2 + phi_in^2)) + 4*phi_inpp; 

    %% Init Jac (non-end points)
    one = ones(m,1);
    Lapl = -0.5 * spdiags([one, -2*one, one],[-1,0,1],m,m) / h2;
    A11 = Lapl + spdiags(3*phi_r.^2+phi_i.^2+comm,0,m,m);
    A12 = 2 * spdiags(phi_r.*phi_i,0,m,m);
    A21 = A12;
    A22 = Lapl + spdiags(3*phi_i.^2+phi_r.^2+comm,0,m,m);

    %% BCs
    A11(1,1) = -0.5 * ( dOmega_l_X0 * phi_r(1) + Omega_l(1) ) + 3 * phi_rr2(1) + phi_ii2(1) + comm(1); %***
    A11(1,2) = -0.5 * dOmega_l_X1 * phi_r(1); 
    A11(1,3) = -0.5 * dOmega_l_X2 * phi_r(1);
    A11(m,m-2) = -0.5 * dOmega_r_XNm1 * phi_r(m);
    A11(m,m-1) = -0.5 * dOmega_r_XN * phi_r(m);
    A11(m,m) = - 0.5 * (dOmega_r_XNp1*phi_r(m) + Omega_r) + 3 * phi_rr2(m) + phi_ii2(m) + comm(m); 
    
    A22(1,1) = -0.5 * ( dOmega_l_Y0 * phi_i(1) + Omega_l ) + 3 * phi_ii2(1) + phi_rr2(1) + comm(1);
    A22(1,2) = -0.5 * ( dOmega_l_Y1 * phi_i(1) );
    A22(1,3) = -0.5 * ( dOmega_l_Y2 * phi_i(1) );
    A22(m,m-2) = -0.5 * dOmega_r_YNm1 * phi_i(m);
    A22(m,m-1) = -0.5 * dOmega_r_YN * phi_i(m);
    A22(m,m) = -0.5 * (dOmega_r_YNp1 * phi_i(m) + Omega_r) + 3 * phi_ii2(m) + phi_rr2(m) + comm(m);

    A12(1,1) = -0.5 * dOmega_l_Y0 * phi_r(1) + 2 * phi_r(1) * phi_i(1);
    A12(1,2) = -0.5 * dOmega_l_Y1 * phi_r(1);
    A12(1,3) = -0.5 * dOmega_l_Y2 * phi_r(1);
    A12(m,m-2) = -0.5 * dOmega_r_YNm1 * phi_r(m);
    A12(m,m-1) = -0.5 * dOmega_r_YN * phi_r(m);
    A12(m,m) = -0.5 * dOmega_r_YNp1*phi_r(m) + 2 * phi_r(m)*phi_i(m);
    
    A21(1,1) = -0.5 * dOmega_l_X0 * phi_i(1) + 2 * phi_r(1) * phi_i(1);
    A21(1,2) = -0.5 * dOmega_l_X1 * phi_i(1);
    A21(1,3) = -0.5 * dOmega_l_X2 * phi_i(1);
    A21(m,m-2) = -0.5 * dOmega_r_XNm1 * phi_i(m);
    A21(m,m-1) = -0.5 * dOmega_r_XN * phi_i(m);
    A21(m,m) = -0.5 * dOmega_r_XNp1 * phi_i(m) + 2 * phi_r(m) * phi_i(m);

    %% Construct jac
    jac = [A11 A12; A21 A22];

end
    
