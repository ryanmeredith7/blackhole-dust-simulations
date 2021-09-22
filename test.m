% Paramaters to set up the initial conditions and solver.
p0 = 0.02;
x0 = 10;
a0 = 10;
dx = 0.01;
dt = 0.001;
n = 30000;

% Midpoints of the cells.
xs = (dx/2:dx:a0).';

% Creates the initial condition vectors, first the outside case, the the inside
% case.
ai = (x0 / a0 ./ xs) .^ 2;
bi = -asin(sqrt(8*pi*p0/3*x0^3 ./ xs .^ 3 - ai));

ai(xs < x0) = a0 ^ -2;
bi(xs < x0) = -asin(sqrt(8*pi*p0/3 - a0 ^ -2));

while true

end
