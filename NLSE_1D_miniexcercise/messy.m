function sol = nls_onecomp_fdm(interval, bv, m, nstep, mu, gam)
    % Finite difference method to solve eq. 1.3 using Newton's method
    % Inputs:  boundary conditions and number of grid points
    % OUtputs: Approximate solution to eq. 1.3

    % Make boundary conditions more ascessible
    a = interval(1);
    b = interval(2);

    % Initialize for Newton's
    h = (b-a)/(m+1);
    if gam == -1 % init guesses
%         bright = @(x) sech(x);
%         sol = bright(linspace(-20,20,m))';
        
        sol = zeros(m,1);
%         sol(199) = 1;  
% 
%         l1 = @(x) 0.2*x+1;
%         l2 = @(x) -0.2*x+1;
%         sol(150:199) = l1(linspace(-5,0,50));
%         sol(200:250) = l2(linspace(0,5,51));
        
%         sol = (0.5*cos(linspace(-20,20,m).*pi./(20)) + 0.5)';

        sol(150:250) = (0.5*cos(linspace(-5,5,101).*pi./(5)) + 0.5)';
        
    %     sol(175:225) = 0.5;
%         sol = ones(m,1);
    elseif gam == 1
%         dark = @(x) tanh(x);
%         sol = dark(linspace(-20,20,m))';
        sol = zeros(m,1);
        sol(1:198) = -1;       
        sol(200:end) = 1;
    end
    
    % Newton's method step
    if gam == -1
        for i = 1:nstep
            sol = sol - (1/h^2)*jacobian(sol,interval, m, mu, gam)\v_bright(sol,interval,bv,m, mu, gam);
%             sol = sol - jacobian(sol,interval, m, mu, gam)\v_bright(sol,interval,bv,m, mu, gam);
        end
        
    elseif gam == 1
        for i = 1:nstep
            sol = sol - (1/h^2)*jacobian(sol,interval, m, mu, gam)\v_dark(sol,interval,bv,m, mu, gam);
%             sol = sol - jacobian(sol,interval, m, mu, gam)\v_dark(sol,interval,bv,m, mu, gam);
        end
    end

end

function J = jacobian(sol,interval, m, mu, gam)
    % Makes Jacobian matrix

    h = (interval(1)-interval(2))/(m+1);
    J = zeros(m,m);
    
    % main diagonal
    for i = 1:m
%         J(i,i) = -2 + 2*mu*h^2 - 6*gam*h^2*(sol(i)^2); %taking deriv
%         J(i,i) = -2 - 6*h^2*gam*sol(i)^2 - mu; %*** kent
%         J(i,i) = -2*(1/h^2) - 6*gam*sol(i)^2 - mu*(1/h^2); %doesnt break with tru input
        J(i,i) = -2*(1/h^2) - 6*gam*sol(i)^2 + 2*mu; % rederived
%         J(i,i) = -2 - 2*(h^2)*gam*sol(i)^2 + 2*mu*(h^2) - 4*(h^2)*gam*sol(i);
%         J(i,i) = 1 + h^2 *3*gam*sol(i)^2 - mu;
%         J(i,i) = -2 - 6*h^2*gam*sol(i)^2 + mu;
    end
    
        
%     J(1,1) = -1;
%     J(1,2) = 0;
%     J(2,1) = 0;
%     
%     J(m,m) = -1;
%     J(m,m-1) = 0;
%     J(m-1,m) = 0;
    
    
    % other diagonals of tridiag
    for i = 1:m-1
        if gam == -1
            J(i,i+1) = 1;
            J(i+1,i) = 1;
        elseif gam == 1
            J(i,i+1) = 1/h^2;
            J(i+1,i) = 1/h^2;
        end
    end
end


function y = v_bright(sol,interval,bv,m, mu, gam)
    % Eq. 1.3 vector

    % Init
    y = zeros(m,1);
    h = (interval(2)-interval(1))/(m+1);

    % Lower bv
%     y(1) = bv(1) + (-2+2*mu*h^2)*sol(1) - 2*gam*h^2*sol(1)^3 + sol(2);
    y(1) = (sol(2) - 2*sol(1) + bv(1))*1/(h^2) + -2*(gam*sol(1)^2 - mu)*sol(1);
%     y(1) = (sol(2) - 2*sol(1) + bv(1)) + -2*(h^2)*(gam*sol(1)^2 - mu)*sol(1);
    
    % Upper bv
%     y(m) = sol(m-1) + (-2+2*mu*h^2)*sol(m) - 2*gam*h^2*sol(m)^3 + bv(2);
    y(m) = (bv(2) - 2*sol(m) + sol(m-1))*1/(h^2) + -2*(gam*sol(m)^2 - mu)*sol(m);
%     y(m) = (bv(2) - 2*sol(m) + sol(m-1)) + -2*(h^2)*(gam*sol(m)^2 - mu)*sol(m);

    % Intermediate bv
    for i=2:m-1
%         y(i) = sol(i-1) + (-2+2*mu*h^2)*sol(i) - 2*gam*h^2*sol(i)^3 + sol(i+1);
        y(i) = (sol(i+1) - 2*sol(i) + sol(i-1))*1/(h^2) + -2*(gam*sol(i)^2 - mu)*sol(i);
%         y(i) = (sol(i+1) - 2*sol(i) + sol(i-1))+ -2*(h^2)*(gam*sol(i)^2 - mu)*sol(i);
    end
end


function y = v_dark(sol,interval,bv,m, mu, gam)
    % Eq. 1.3 vector

    % Init
    y = zeros(m,1);
    h = (interval(2)-interval(1))/(m+1);

    % Lower bv
%     y(1) = bv(1) + (-2+2*mu*h^2)*sol(1) - 2*gam*h^2*sol(1)^3 + sol(2);
    y(1) = (sol(2) - 2*sol(1) + sol(2))*1/(h^2) + -2*(gam*sol(1)^2 - mu)*sol(1);
%     y(1) = (sol(2) - 2*sol(1) + bv(1)) + -2*(h^2)*(gam*sol(1)^2 - mu)*sol(1);
    
    % Upper bv
%     y(m) = sol(m-1) + (-2+2*mu*h^2)*sol(m) - 2*gam*h^2*sol(m)^3 + bv(2);
    y(m) = (sol(m-1) - 2*sol(m) + sol(m-1))*1/(h^2) + -2*(gam*sol(m)^2 - mu)*sol(m);
%     y(m) = (bv(2) - 2*sol(m) + sol(m-1)) + -2*(h^2)*(gam*sol(m)^2 - mu)*sol(m);

    % Intermediate bv
    for i=2:m-1
%         y(i) = sol(i-1) + (-2+2*mu*h^2)*sol(i) - 2*gam*h^2*sol(i)^3 + sol(i+1);
        y(i) = (sol(i+1) - 2*sol(i) + sol(i-1))*1/(h^2) + -2*(gam*sol(i)^2 - mu)*sol(i);
%         y(i) = (sol(i+1) - 2*sol(i) + sol(i-1))+ -2*(h^2)*(gam*sol(i)^2 - mu)*sol(i);
    end
end