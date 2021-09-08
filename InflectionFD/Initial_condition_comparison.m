%% Code used to compare errors between orders of approximation of initial conditions

%Default parameter definitions
xmax = 5;
t0 = -10;
tmax = 10;
dx = 0.001;
dt = 0.001;
WG = 1;
ICorder = 2;
t0vec = linspace(-5,-1)';
x=(0:dx:xmax);
M=length(x)-1;
t=(t0:dt:tmax);
N=length(t)-1;
nu = -airy0(0, WG);
c1 = -1i * nu^2 * 32/225;
c2 = 32/45 * nu - 512/50625 * nu^4;
X = x*(-2*t0).^(1/3);
phi0 = airy(X-nu);
phi1 = (c1 + 1i/3 .* X.^2) .* airy(X-nu);
phi2 = (c2 -16/45 .* X + 64/1350 * nu^2 * X.^2 -1/18 * X.^4) .* airy(X-nu) + (64/135 * nu * X + 16/45 * X.^2) .* airy(1,(X-nu));
u0=(-t0)^(1/6)*exp(1i*nu*3/5*2^(-1/3)*(-t0)^(5/3))*(phi0);
u1=(-t0)^(1/6)*exp(1i*nu*3/5*2^(-1/3)*(-t0)^(5/3))*(phi0 + (-2 * t0)^(-5/3) * phi1);
u2=(-t0)^(1/6)*exp(1i*nu*3/5*2^(-1/3)*(-t0)^(5/3))*(phi0 + (-2 * t0)^(-5/3) * phi1 + ((-2 * t0)^(-10/3)) * phi2);

% figure;
% plot1r = plot(x,real(u0), 'k-', 'DisplayName', 'Real Part'); hold on;
% plot1i = plot(x,imag(u0), 'r-', 'DisplayName', 'Imaginary Part');
% xlabel('X');
% ylabel('\psi(X,t_0)');
% legend;
% title('Initial Condition on t=t_0 up to Zeroth Order');
% 
% figure;
% plot2r = plot(x,real(u1), 'k-', 'DisplayName', 'Real Part'); hold on;
% plot2i = plot(x,imag(u1), 'r-', 'DisplayName', 'Imaginary Part');
% xlabel('X');
% ylabel('\psi(X,t_0)');
% legend;
% title('Initial Condition on t=t_0 up to First Order');
% 
% figure;
% plot3r = plot(x,real(u2), 'k-', 'DisplayName', 'Real Part'); hold on;
% plot3i = plot(x,imag(u2), 'r-', 'DisplayName', 'Imaginary Part');
% xlabel('X');
% ylabel('\psi(X,t_0)');
% legend;
% title('Initial Condition on t=t_0 up to Second Order');


% figure;
% plot1r = plot(x,(real(u1)-real(u2)).^2, 'k-', 'DisplayName', 'Real Part'); hold on;
% plot1i = plot(x,(imag(u1)-imag(u2)).^2, 'r-', 'DisplayName', 'Imaginary Part');
% xlabel('X');
% ylabel('Square error in initial conditions at t_0');
% legend;
% title('Error between first and second order ICs, t_0 = -10');
% 
% figure;
% plot1r = plot(x,(real(u1)-real(u0)).^2, 'k-', 'DisplayName', 'Real Part'); hold on;
% plot1i = plot(x,(imag(u1)-imag(u0)).^2, 'r-', 'DisplayName', 'Imaginary Part');
% xlabel('X');
% ylabel('Square error in initial conditions at t_0');
% legend;
% title('Error between zeroth and first order ICs, t_0 = -10');

u0mat=(-t0vec).^(1/6).*exp(1i*nu*3/5*2^(-1/3).*(-t0vec).^(5/3)).*(phi0);
u1mat=(-t0vec).^(1/6).*exp(1i*nu*3/5*2^(-1/3).*(-t0vec).^(5/3)).*(phi0 + (-2 * t0vec).^(-5/3) * phi1);
u2mat=(-t0vec).^(1/6).*exp(1i*nu*3/5*2^(-1/3).*(-t0vec).^(5/3)).*(phi0 + (-2 * t0vec).^(-5/3) * phi1 + ((-2 * t0vec).^(-10/3)) * phi2);

zfdiffreal = real((real(u0mat)-real(u1mat)).^2);
fsdiffreal = real((real(u1mat)-real(u2mat)).^2);
zfdiffimag = (imag(u0mat)-imag(u1mat)).^2;
fsdiffimag = (imag(u1mat)-imag(u2mat)).^2;

[t0t0,xx] = meshgrid(t0vec, x);

figure;
surfrealzfdiff = surf(t0t0,xx,real(zfdiffreal')); hold on;
surfrealzfdiff.EdgeColor = 'none';
ylim([0,5]);
xlabel('t0');
ylabel('X');
zlabel('Square error in initial conditions at t_0');
title('Error between real part of zeroth and first order ICs');

figure;
surfimagzfdiff = surf(t0t0,xx,real(zfdiffimag')); hold on;
surfimagzfdiff.EdgeColor = 'none';
ylim([0,5]);
xlabel('t0');
ylabel('X');
zlabel('Square error in initial conditions at t_0');
title('Error between imaginary part of zeroth and first order ICs');

figure;
surfrealfsdiff = surf(t0t0,xx,real(fsdiffreal')); hold on;
surfrealfsdiff.EdgeColor = 'none';
ylim([0,5]);
xlabel('t0');
ylabel('X');
zlabel('Square error in initial conditions at t_0');
title('Error between real part of first and second order ICs');

figure;
surfimagfsdiff = surf(t0t0,xx,real(fsdiffimag')); hold on;
surfimagfsdiff.EdgeColor = 'none';
ylim([0,5]);
xlabel('t0');
ylabel('X');
zlabel('Square error in initial conditions at t_0');
title('Error between imaginary part of first and second order ICs');


