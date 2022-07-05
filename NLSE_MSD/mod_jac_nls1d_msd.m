function J = mod_jac_nls1d_msd( psi, params )

% Unpack parameters:
npts = params.m;
   dx = params.h;
    V = params.Vstat;
 mu = params.mu;

% Decompose the field into real and imaginary parts + auxiliary calculations:
X = psi(1:npts); X2 = X.^2; 
Y = psi(npts+1:2*npts); Y2 = Y.^2;

% Auxiliary computations:
dx2 = dx * dx;
XX0 = X(1); XX1 = X(2); XX2 = X(3);
YY0 = Y(1); YY1 = Y(2); YY2 = Y(3);
  V1 = V(2); V0 = V(1);

XXNp1 = X(npts); XXN = X(npts-1); XXNm1 = X(npts-2); 
YYNp1 = Y(npts); YYN = Y(npts-1); YYNm1 = Y(npts-2);
      VN = V(npts-1); VNp1 = V(npts);

     comm =  V - mu*ones(npts,1); 
    % term_l = ((XX2 - 2*XX1 + XX0)*XX1 + (YY2 - 2*YY1 + YY0)*YY1)...
    %             /(dx2*(XX1^2 + YY1^2)) - 2*(XX1^2 + YY1^2 - XX0^2 - YY0^2 + V1 - V0);
    % term_r = ((XXNp1 - 2*XXN + XXNm1)*XXN + (YYNp1 - 2*YYN + YYNm1)*YYN) ...
     %           /(dx2*(XXN^2 + YYN^2)) - 2*(XXN^2 + YYN^2 - XXNp1^2 - YYNp1^2 + VN - VNp1);
  %   term_l = dx2.^(-1).*(X1.^2+Y1.^2).^(-1).*(X1.*(X0+(-2).*X1+X2)+Y1.*(Y0+(-2).*Y1+Y2));
%     term_l = ( ( X(3)-2*X(2)+X(1) ) * X(2) + ( Y(3)-2*Y(2)+Y(1) ) * Y(2) ) / ( dx^2 * ( X(2)^2+Y(2)^2 ) );
   %  term_r = ( ( X(npts)-2*X(npts-1)+X(npts-2) ) * X(npts-1) + ( Y(npts)-2*Y(npts-1)+Y(npts-2) ) * Y(npts-1) )...
    %/ ( dx^2 * ( X(npts-1)^2+Y(npts-1)^2 ) );

% Omega_l and its derivatives:
%        Omega_l = term_l - 2 * ( X(2).^2 + Y(2).^2 - X(1).^2 - Y(1).^2 + V(2) - V(1) );
        Omega_l = ((XX2 - 2*XX1 + XX0)*XX1 + (YY2 - 2*YY1 + YY0)*YY1)...
                 /(dx2*(XX1^2 + YY1^2)) - 2*(XX1^2 + YY1^2 - XX0^2 - YY0^2 + V1 - V0);
dOmega_l_X0 = 4*XX0 + XX1/(dx2*(XX1^2 + YY1^2)); %4 * X(1) + X(2) / ( dx^2 * ( X(2).^2 + Y(2).^2 ) );
dOmega_l_X1 = -4*XX1 + (XX0 - 4*XX1 + XX2)/(dx2*(XX1^2 + YY1^2)) ...
                         - ( 2*XX1*(XX1*(XX0 - 2*XX1 + XX2) + YY1*(YY0 - 2*YY1 + YY2)))/(dx2*(XX1^2 + YY1^2)^2);
                           %-4 * X(2) + ( X(1) - 4 * X(2) + X(3) ) / ( dx^2 * ( X(2).^2 + Y(2).^2 ) ) ...
                           %-2 * X(2) * ( X(2) * ( X(1)-2*X(2)+X(3) ) + Y(2) * ( Y(1)-2*Y(2)+Y(3) ) ) ...
                           %/ ( dx^2 * ( X(2).^2 + Y(2).^2 ).^2 );
dOmega_l_X2 = XX1/(dx2*(XX1^2 + YY1^2)); %X(2) / ( dx^2 * ( X(2).^2 + Y(2).^2 ) );
dOmega_l_Y0 = 4*YY0 + YY1/(dx2*(XX1^2 + YY1^2)); %4 * Y(1) + Y(2) / ( dx^2*( X(2).^2+Y(2).^2 ) );
dOmega_l_Y1 = -4*YY1 + (YY0 - 4*YY1 + YY2)/(dx2*(XX1^2 + YY1^2)) ...
                           - (2*YY1*(XX1*(XX0 - 2*XX1 + XX2) + YY1*(YY0 - 2*YY1 + YY2)))...
                            /(dx2*(XX1^2 + YY1^2)^2); % -4 * Y(2) + ( Y(1) - 4 * Y(2) + Y(3) )/(dx^2*( X(2).^2+Y(2).^2) ) ...
                          %- 2 * Y(2) * ( X(2) * ( X(1)-2*X(2)+X(3) ) + Y(2) * ( Y(1) - 2*Y(2) + Y(3) ) )...
                          %/( dx^2 * ( X(2).^2+Y(2).^2 ).^2 ); 
dOmega_l_Y2 = YY1/(dx2*(XX1^2 + YY1^2)); %Y(2) / ( dx^2 * ( X(2).^2 + Y(2).^2 ) );

% Omega_r and its derivatives:
%              Omega_r = term_r - 2 * ( X(npts-1).^2 + Y(npts-1).^2 - X(npts).^2 - Y(npts).^2 + V(npts-1) - V(npts) );
              Omega_r = ((XXNp1 - 2*XXN + XXNm1)*XXN + (YYNp1 - 2*YYN + YYNm1)*YYN) ...
                     /(dx2*(XXN^2 + YYN^2)) - 2*(XXN^2 + YYN^2 - XXNp1^2 - YYNp1^2 + VN - VNp1);
dOmega_r_XNm1 = XXN/(dx2*(XXN^2 + YYN^2)); %X(npts-1) / ( dx^2 * ( X(npts-1).^2 + Y(npts-1).^2 ) );
     dOmega_r_XN = -4*XXN + (-4*XXN + XXNm1 + XXNp1)/(dx2*(XXN^2 + YYN^2))...
                                 - (2*XXN*(XXN*(-2*XXN + XXNm1 + XXNp1) ...
                                 + YYN * (-2*YYN + YYNm1 + YYNp1)))...
                                 /(dx2*(XXN^2 + YYN^2)^2);
                               %-4 * X(npts-1) + ( -4*X(npts-1)+X(npts-2) + X(npts) ) ...
                               %/ ( dx^2*( X(npts-1).^2+Y(npts-1).^2 ) )...
                              %-2*X(npts-1) * ( X(npts-1)*( -2*X(npts-1) + X(npts-2) + X(npts) ) + ...
                              %Y(npts-1) * ( -2 * Y(npts-1) + Y(npts-2) + Y(npts) ) ) ...
                              %/ (dx^2 * ( X(npts-1).^2 + Y(npts-1).^2 ).^2);
 dOmega_r_XNp1 = 4*XXNp1 + XXN/(dx2*(XXN^2 + YYN^2)); %4 * X(npts) + X(npts-1) / ( dx^2 * ( X(npts-1).^2 + Y(npts-1).^2 ) );
dOmega_r_YNm1 = YYN/(dx2*(XXN^2 + YYN^2)); %Y(npts-1) / ( dx^2 * ( X(npts-1).^2 + Y(npts-1).^2 ) );
     dOmega_r_YN = -4*YYN + (-4*YYN + YYNm1 + YYNp1)/(dx2*(XXN^2 + YYN^2)) ...
                                 - (2*YYN*(XXN*(-2*XXN + XXNm1 + XXNp1) ...
                                 + YYN*(-2*YYN + YYNm1 + YYNp1)))/(dx2*(XXN^2 + YYN^2)^2);
                                %-4*Y(npts-1)+(-4*Y(npts-1)+Y(npts-2)+Y(npts))/(dx^2*(X(npts-1).^2+Y(npts-1).^2))...
                                % -2*Y(npts-1) * ( X(npts-1)*(-2*X(npts-1) + X(npts-2) +X(npts) ) ...
                                 %+Y(npts-1)*(-2*Y(npts-1)+Y(npts-2)+Y(npts)) ) ...
                                 %/ (dx^2*( X(npts-1).^2 + Y(npts-1).^2 ).^2);
dOmega_r_YNp1 = YYN/(dx2*(XXN^2 + YYN^2)) + 4*YYNp1; %4 * Y(npts) + Y(npts-1) / ( dx^2*( X(npts-1).^2+Y(npts-1).^2) );

% Pre-allocate stuff + preliminary computations:
 one = ones(npts,1);
Lapl = -0.5 * spdiags([one, -2*one, one],[-1,0,1],npts,npts) / dx^2; % The Laplacian with 0 BCs.
 J11 = Lapl + spdiags(3*X.^2+Y.^2+comm,0,npts,npts);
 J12 = 2 * spdiags(X.*Y,0,npts,npts);
 J21 = J12;
 J22 = Lapl + spdiags(3*Y.^2+X.^2+comm,0,npts,npts);

% J11 block:
% -------------------- left boundary
J11(1,1) = -0.5 * ( dOmega_l_X0 * X(1) + Omega_l(1) ) + 3 * X2(1) + Y2(1) + comm(1);
J11(1,2) = -0.5 * dOmega_l_X1 * X(1); 
J11(1,3) = -0.5 * dOmega_l_X2 * X(1);
% -------------------- right boundary
J11(npts,npts-2) = -0.5 * dOmega_r_XNm1 * X(npts);
J11(npts,npts-1) = -0.5 * dOmega_r_XN * X(npts);
J11(npts,npts) = - 0.5 * (dOmega_r_XNp1*X(npts) + Omega_r) + 3 * X2(npts) + Y2(npts) + comm(npts); 

% J12 block:
% -------------------- left boundary
J12(1,1) = -0.5 * dOmega_l_Y0 * X(1) + 2 * X(1) * Y(1);
J12(1,2) = -0.5 * dOmega_l_Y1 * X(1);
J12(1,3) = -0.5 * dOmega_l_Y2 * X(1);
% -------------------- right boundary
J12(npts,npts-2) = -0.5 * dOmega_r_YNm1 * X(npts);
J12(npts,npts-1) = -0.5 * dOmega_r_YN * X(npts);
   J12(npts,npts) = -0.5 * dOmega_r_YNp1*X(npts) + 2 * X(npts)*Y(npts);

% J21 block:
% -------------------- left boundary
J21(1,1) = -0.5 * dOmega_l_X0 * Y(1) + 2 * X(1) * Y(1);
J21(1,2) = -0.5 * dOmega_l_X1 * Y(1);
J21(1,3) = -0.5 * dOmega_l_X2 * Y(1);
% -------------------- right boundary
J21(npts,npts-2) = -0.5 * dOmega_r_XNm1 * Y(npts);
J21(npts,npts-1) = -0.5 * dOmega_r_XN * Y(npts);
J21(npts,npts) = -0.5 * dOmega_r_XNp1 * Y(npts) + 2 * X(npts) * Y(npts);

% J22 block:
% -------------------- left boundary
J22(1,1) = -0.5 * ( dOmega_l_Y0 * Y(1) + Omega_l ) + 3 * Y2(1) + X2(1) + comm(1);
J22(1,2) = -0.5 * ( dOmega_l_Y1 * Y(1) );
J22(1,3) = -0.5 * ( dOmega_l_Y2 * Y(1) );
% -------------------- right boundary
J22(npts,npts-2) = -0.5 * dOmega_r_YNm1 * Y(npts);
J22(npts,npts-1) = -0.5 * dOmega_r_YN * Y(npts);
   J22(npts,npts) = -0.5 * (dOmega_r_YNp1 * Y(npts) + Omega_r) + 3 * Y2(npts) + X2(npts) + comm(npts);

% Return the "PRAGMA":
J = [J11, J12; ... 
       J21, J22 ];

end