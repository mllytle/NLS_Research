function jac = jac_32(phi, params)
% x is concatenation of R and Imag comps of wave eq. at disc. points
    
    %% get params
    s = params.s;
    m = params.m;
    a = params.a;
    Vfunct = params.V;
    V = Vfunct(phi);
    mu = params.mu;
    h = params.h;

    phi_r = phi(1:m);
    phi_r2 = phi_r.^2;
    phi_i = phi(m+1:2*m);
    phi_i2 = phi_i.^2;

    A11 = zeros(m);
    A12 = zeros(m);
    A21 = zeros(m);
    A22 = zeros(m);

    % construct leading and trailing rows
    %% R eq w/t respect to R
    
    % first row bc
    %   off diag

%     dAdphiR_1 = phi_r(2)/(h^2*(phi_r(2)^2 + phi_i(2)^2));
%     A_R1 = 
% 
%     A11(1,1) = phi_r(1)*(dAdphiR_1 + (2*s/a)*phi_r(1)) + ...
%         A_R1 + (1/a)*(s*(phi_r(2)^2 + phi_i(2)^2) - V(2) ...
%         - s*(phi_r(1)^2 + phi_i(1)^2) + V(1));

%     A11(1,3) = phi_r(3)*(dAdphiR_1 + (2*s/a)*phi_r(3)) + ...
%         A + (1/a)*(s*(phi_r(2)^2 + phi_i(2)^2) - V(2) ...
%         - s*(phi_r(3)^2 + phi_i(3)^2) + V(3));

%     A11(1,2) = 

    A11(1,1) = -s*(3*phi_r(1)^2+phi_i(1)^2) + (V(1)-mu) ... % *** s + or -?
        -a*((2*phi_r(1)*phi_r(2))/(h^2*(phi_r(2)^2 + phi_i(2)^2)) + (1/a)*(s*phi_r(2)^2 ...
        -s*3*phi_r(1)^2) + V(1) - V(2)); 
    A11(1,2) = -a*(-2*phi_r(1)/(h^2))*(2*phi_r(2)*phi_i(2)^2/(phi_r(2)^2 + ...
        phi_i(2)^2)^2);
    A11(1,3) = -a*(1/(h^2*(phi_r(2)^2 + phi_i(2)^2)))*phi_r(1);

    
    A11(m,m) = -s*(3*phi_r(m)^2+phi_i(m)^2) + (V(m)-mu) ... 
        -a*((2*phi_r(m)*phi_r(m-1))/(h^2*(phi_r(m-1)^2 + phi_i(m-1)^2)) + (1/a)*(s*phi_r(m-1)^2-s*3*phi_r(m)^2)...
        + V(m) - V(m-1)); 
    A11(m,m-1) = -a*(-2*phi_r(m)/(h^2))*(2*phi_r(m-1)*phi_i(m-1)^2/(phi_r(m-1)^2 + phi_i(m-1)^2)^2);
    A11(m,m-2) = -a*(1/(h^2*(phi_r(m-1)^2 + phi_i(m-1)^2)))*phi_r(m);
    
    %% I eq w/t respect to I
    A22(1,1) = -s*(3*phi_i(1)^2+phi_r(1)^2) + (V(1)-mu) ...
        -a*((2*phi_i(1)*phi_i(2))/(h^2*(phi_r(2)^2 + phi_i(2)^2)) + ...
        (1/a)*(s*phi_i(2)^2-s*3*phi_i(1)^2) + V(1) - V(2)); 
    A22(1,2) = -a*(-2*phi_i(1)/(h^2))*(2*phi_i(2)*phi_r(2)^2/(phi_i(2)^2 + phi_r(2)^2)^2);
    A22(1,3) = -a*(1/(h^2*(phi_i(2)^2 + phi_r(2)^2)))*phi_i(1);

    A22(m,m) = -s*(3*phi_i(m)^2+phi_r(m)^2) + (V(m)-mu) ... 
        -a*((2*phi_i(m)*phi_i(m-1))/(h^2*(phi_r(m-1)^2 + phi_i(m-1)^2)) +...
        (1/a)*(s*phi_i(m-1)^2-s*3*phi_i(m)^2) + V(m) - V(m-1));
    A22(m,m-1) = -a*(-2*phi_i(m)/(h^2))*(2*phi_i(m-1)*phi_r(m-1)^2/(phi_i(m-1)^2 + phi_r(m-1)^2)^2);
    A22(m,m-2) = -a*(1/(h^2*(phi_i(m-1)^2 + phi_r(m-1)^2)))*phi_i(m);
    
    %% R eq w/t respect to I
    A12(1,1) = -s*2*phi_i(1)*phi_r(1) -a*(phi_i(2)/(h^2*(phi_r(2)^2 + ...
        phi_i(2)^2)))*phi_r(1);
    A12(1,2) = -s*(-4*phi_i(2)/(h^2*(phi_r(2)^2+phi_i(2)^2)))*phi_r(1);
    A12(1,3) = -s*(phi_i(2)/(h^2*(phi_r(2)^2+phi_i(2)^2)))*phi_r(1);

    A12(m,m) = -s*2*phi_i(m)*phi_r(m) -a*(phi_i(m-1)/(h^2*(phi_r(m-1)^2 + ...
        phi_i(m-1)^2)))*phi_r(m);
    A12(m,m-1) = -s*(-4*phi_i(m-1)/(h^2*(phi_r(m-1)^2+phi_i(m-1)^2)))*phi_r(m);
    A12(m,m-2) = -s*(phi_i(m-1)/(h^2*(phi_r(m-1)^2+phi_i(m-1)^2)))*phi_r(m);

    %% I eq w/t respect to R 
    A21(1,1) = -s*2*phi_i(1)*phi_r(1) -a*(phi_r(2)/(h^2*(phi_r(2)^2 + ...
        phi_i(2)^2)))*phi_i(1);
    A21(1,2) = -s*(-4*phi_r(2)/(h^2*(phi_r(2)^2+phi_i(2)^2)))*phi_i(1);
    A21(1,3) = -s*(phi_r(2)/(h^2*(phi_r(2)^2+phi_i(2)^2)))*phi_i(1);

    A21(m,m) = -s*2*phi_i(m)*phi_r(m) -a*(phi_r(m-1)/(h^2*(phi_r(m-1)^2 + ...
        phi_i(m-1)^2)))*phi_i(m);
    A21(m,m-1) = -s*(-4*phi_r(m-1)/(h^2*(phi_r(m-1)^2+phi_i(m-1)^2)))*phi_i(m);
    A21(m,m-2) = -s*(phi_r(m-1)/(h^2*(phi_r(m-1)^2+phi_i(m-1)^2)))*phi_i(m);


    %% construct diagonal
    for ii = 2:m-1
        A11(ii,ii) = (-2*a/h^2) - s*(3*phi_r(ii)^2+phi_i(ii)^2) + (V(ii)-mu);
        A12(ii,ii) = -s*2*phi_i(ii)*phi_r(ii);
        A22(ii,ii) = (-2*a/h^2) - s*(phi_r(ii)^2+3*phi_i(ii)^2) + (V(ii)-mu);
        A21(ii,ii) = -s*2*phi_i(ii)*phi_r(ii);
    end
    
    %% construct off diagonals
    for kk = 1:m-1
        A11(kk,kk+1) = -a/h^2;
        A11(kk+1,kk) = -a/h^2;

        A22(kk,kk+1) = -a/h^2;
        A22(kk+1,kk) = -a/h^2;

        A12(kk,kk+1) = 0;
        A12(kk+1,kk) = 0;

        A21(kk,kk+1) = 0;
        A21(kk+1,kk) = 0;
    end

    % each quadrant N+2 x N+2    

    %% Construct jac
    jac = [A11 A12; A21 A22];

end
    