% Paramaters to set up the initial conditions and solver.
p0 = 0.02;
x0 = 4;
a0 = 5;
dx = 0.01;
dt = 0.001;
n = uint64(30000);

% Midpoints of the cells.
xs = (dx/2:dx:a0).';

% Creates the initial condition vectors, first the outside case, then the
% inside case.
ai = (x0 / a0 ./ xs) .^ 2;
bi = -asin(sqrt(8*pi*p0/3*x0^3 ./ xs .^ 3 - ai));

ai(xs < x0) = a0 ^ -2;
bi(xs < x0) = -asin(sqrt(8*pi*p0/3 - a0 ^ -2));

% Calls the function that solves the equations.
[a, b] = solver(ai, bi, 0, dx, dt, n);

% Calculates the denssity.
bmid = (b(:,1:end-1) + b(:,2:end)) ./ 2;
p = (diff(b, 1, 2) + diff(a, 1, 2) ./ sin(2 .* bmid)) ./ (-4*pi*dt);

speed = 20;

% Plays a short movie of one of the solution values.
figure(Name="Animation of a");
animate(xs, a, speed);
uiwait(msgbox("Press OK to continue.", "Done plotting a"));

figure(Name="Animation of alpha");
animate(xs, 1 - xs .^ 2 .* a, speed);
uiwait(msgbox("Press OK to continue.", "Done plotting alpha"));

figure(Name="Animation of beta");
animate(xs, b, speed);
uiwait(msgbox("Press OK to continue.", "Done plotting beta"));

figure(Name="Animation of rho");
animate(xs, p, speed, [0, 0.2]);
uiwait(msgbox("Press OK to continue.", "Done plotting rho"));
