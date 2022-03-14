% Paramaters to set up the initial conditions and solver.
p0 = 1e-5;
x0 = 40;
a0 = 50;
dx = 0.1;
dt = 0.02;

% Cell boundries
x = (0:dx:a0).';
xo2 = x ./ 2;
% Cell midpoints
xm = x(2:end) - dx/2;

% Creates the initial condition vectors, first the outside case, then the inside
% case.
b = -asin(sqrt(8*pi*p0/3*x0^3 ./ xm .^ 3));
b(xm < x0) = -asin(sqrt(8*pi*p0/3));
a = zeros(size(b));

[t1,t2] = time(a, b, 0, dx, dt);
disp([t1,t2]);
