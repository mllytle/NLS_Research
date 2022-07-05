clear
% close all
clc

%% Params
plots = 0;

A = 1;
mu = -1/2;
gam = -1;
m = 399;
h = (40)/(m+1);
% syms u1
% u = u1*ones(m,1);

xdom = linspace(-20,20,m);
u_b = @(x) sech(x);
v = u_b(xdom)';

u_d = @(x) sqrt(mu)*tanh(sqrt(mu)*x);
v2 = u_d(xdom)';

%% Eig
A12 = full(spdiags(gam.*v.^2,0,m,m));
A11 = full(spdiags([-1/(2*h^2)*ones(m,1) gam.*v.*conj(v) -1/(2*h^2)*ones(m,1)],-1:1,m,m));

A = [A11 A12; -conj(A12) -A11];
eigs = eig(A).*(-1i);

%% Plot Eigenvalues for Stability Analysis
figure(1)
plot(eigs, '.')
grid on
xlabel('Real')
ylabel('Imaginary')

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
end



%% Call funct
% w1 = eig_1dnls_b([-20 20], m, 500, mu, gam);
% xdom = linspace(-20,20,length(w1));
% 
% figure(1)
% plot(xdom,w1)

%% Functions
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
    elseif  strcmp(p_or_n, 'p')
        for ll = 1:m
           A(ll,ll) = 1/(2*h^2) + (-2*gam - 1i)*w(ll);
        end
    end
    
    for kk = 1:m-1
        A(kk,kk+1) = 1/(2*h^2);
        A(kk+1,kk) = 1/(2*h^2);
    end
end

%% Appendix

%     l1 = @(x) 0.2*x+1;
%     l2 = @(x) -0.2*x+1;
%     w(150:199) = l1(linspace(-5,0,50));
%     w(200:250) = l2(linspace(0,5,51));
%     w(150:250) = (0.5*cos(linspace(-5,5,101).*pi./(5)) + 0.5)';
%     bright = @(x) sech(1.2*x) + 0.2;
%     bright = @(x) sech(x);
%     w = bright(linspace(-20,20,m))';
%     w2 = bright(linspace(-20,20,m))';
%     dbright = @(x) -2.*(gam.*sech(x).^2 - mu).*sech(x);
%     db = dbright(linspace(-20,20,m));

% A11 = -0.5*dxx + 2*gam*abs(v)^2-mu;
% A12 = gam*v.^2;