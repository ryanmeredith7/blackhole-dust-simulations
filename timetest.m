% Paramaters to set up the initial conditions and solver.
p0 = 0.01;
x0 = 3;
a0 = 5;
dx = 0.01;
dt = 0.001;
n = 30000;

% Cell boundries
x = (0:dx:a0).';
xo2 = x ./ 2;
% Cell midpoints
xm = x(2:end) - dx/2;

% Creates the initial condition vectors, first the outside case, the the inside
% case.
a = (x0 / a0 ./ xm) .^ 2;
b = -asin(sqrt(8*pi*p0/3*x0^3 ./ xm .^ 3 - a));

a(xm < x0) = a0 ^ -2;
b(xm < x0) = -asin(sqrt(8*pi*p0/3 - a0 ^ -2));

[t1,t2] = time(a, b, x0, dx, dt);
disp([t1,t2]);
