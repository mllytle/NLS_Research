clear
close all
clc

%% Params
plots = 1;
plot_eig = 0;
nm_sol = 1;

A = 1;
mu = -1/2;
gam = -1;
m = 399;
h = (40)/(m+1);
% syms u1
% u = u1*ones(m,1);
xdom = linspace(-20,20,m);

params.A = 1;
params.mu = -1/2;
params.gam = -1;
params.m = 399;
params.h = (40)/(m+1);

%% Get sol
if nm_sol
    v = nls_1cbs_fdm([-20 20], m, 500, mu, gam);
    xdom = linspace(-20,20,length(v));
else
    u_b = @(x) sech(x);
    v = u_b(xdom)';
end

%% Eig
A12 = full(spdiags(gam.*v.^2,0,m,m));
A11 = full(spdiags([-1/(2*h^2)*ones(m,1) gam.*v.*conj(v) -1/(2*h^2)*ones(m,1)],-1:1,m,m));

A_eig = [A11 A12; -conj(A12) -A11];
eigs = eig(A_eig)./1i;

%% Plot Eigenvalues for Stability Analysis
if plot_eig
    figure(1)
    plot(eigs, '.')
    grid on
    xlabel('Real')
    ylabel('Imaginary')
end

%% Time Series Plots
if plots
    t = linspace(0,30,m);
    u =@(x,t) exp(-1i*mu.*t).*u_b(x);
    
    for ii = 1:m
        umat(ii,:) = u(xdom,t(ii));
    end
    
    % figure(2)
    % % mesh(xdom,t,umat)
    % mesh(xdom,t,u(xdom,t))
    % grid on
    % xlabel('x')
    % ylabel('t')
    % zlabel('u')
    figure(2)
    hold on
    for ii = 1:m
        clf
        figure(2)
        plot(xdom,abs(u(xdom,t(ii))).^2)
        ylim([-1 1])
        drawnow
    end
    xlabel('x')
    ylabel('u')
    hold off
    grid on
    
    % figure(4)
    % for ii = 1:m
    %     clf
    %     plot3(xdom,real(u(xdom,t(ii))), imag(real(u(xdom,t(ii)))))
    %     ylim([-1 1])
    %     drawnow
    % end
    % xlabel('x')
    % ylabel('u')
    % grid on
    
    % bright = @(x) sech(x);
    % figure()
    % plot(xdom,w1)
    % hold on
    % plot(linspace(-20,20,m),bright(linspace(-20,20,m)))
    % legend('FDM Solution', 'True Solution')

end

%% Functions
function w = nls_1cbs_fdm(inter, m, ns, mu, gam)
    h = (inter(2) - inter(1))/(m+1);
    bright = @(x) sech(1.2*x) + 0.2;
%     bright = @(x) sech(x);
    w = bright(linspace(-20,20,m))';
    for ii = 1:ns
        if ii == 1
%             disp(norm(f(w, m, h, mu, gam)-db')/m)
        end
        w = w - jac(w, m, h, mu, gam)\f(w, m, h, mu, gam); %*** 0.5
%         residual(ii) = norm(f(w, m, h, mu, gam));%norm(w2-w);
    end
end


function A = jac(w, m, h, mu, gam)
    A = zeros(m,m);
    for ll = 1:m
%        A(ll,ll) = (-2 + 2*mu*h^2) - gam*6*(h^2)*w(ll)^2;
       A(ll,ll) = (-2/(h^2) + 2*mu) - gam*6*w(ll)^2;
    end
    
    for kk = 1:m-1
%         A(kk,kk+1)=1;
%         A(kk+1,kk)=1;
        A(kk,kk+1) = 1/(h^2);
        A(kk+1,kk) = 1/(h^2);
    end
end

function w = eig_1dnls_b(inter, m, ns, mu, gam)
    h = (inter(2) - inter(1))/(m+1);
    
    w = zeros(m,1);

    for ii = 1:ns
        if ii == 1
%             disp(norm(f(w, m, h, mu, gam)-db')/m)
        end
        w = w - jac_eig(w, m, h, mu, gam,'p')\f(w, m, h, mu, gam);
%         residual(ii) = norm(f(w, m, h, mu, gam));
    end
    
%     figure(2)
%     plot([1:ns],residual)
%     disp(residual(1))
end

function y = f(w, m, h, mu, gam)
    y = zeros(m,1);
    
    y(1) = w(2) + (-2 - 2*gam*(h^2)*w(1)^2 + 2*mu*h^2)*w(1);
    y(m) = w(m-1) + (-2 - 2*gam*(h^2)*w(m)^2 + 2*mu*h^2)*w(m); 
    
    for jj = 2:m-1
       y(jj) =  w(jj-1) + w(jj+1) + (-2 - 2*gam*(h^2)*w(jj)^2 + 2*mu*h^2)*w(jj);
    end
end

function A = jac_eig(w, m, h, mu, gam, p_or_n)
    A = zeros(m,m);
    
    if strcmp(p_or_n, 'p')
        for ll = 1:m
           A(ll,ll) = 1/(2*h^2) + (-2*gam + 1i)*w(ll);
        end
    elseif  strcmp(p_or_n, 'n')
        for ll = 1:m
           A(ll,ll) = 1/(2*h^2) + (-2*gam - 1i)*w(ll);
        end
    end
    
    for kk = 1:m-1
        A(kk,kk+1) = 1/(2*h^2);
        A(kk+1,kk) = 1/(2*h^2);
    end
end