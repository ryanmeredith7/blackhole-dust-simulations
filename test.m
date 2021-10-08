% Paramaters to set up the initial conditions and solver.
p0 = 0.02;
x0 = 4;
a0 = 5;
dx = 0.01;
dt = 0.001;
n = 30000;

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

% Plays a short movie of one of the solution values.
disp("Showing a");
animate(xs, a);
disp("Done a");
uiwait(gcf);

disp("Showing alpha");
animate(xs, sqrt(1 - xs .^ 2 .* a));
disp("Done alpha");
uiwait(gcf);

disp("Showing beta");
animate(xs, b);
disp("Done beta");
uiwait(gcf);

disp("Showing rho");
animate(xs, p, [0, 0.2]);
disp("Done rho");

function animate(xs, y, l)

    p = plot(xs, y(:,1));
    xlim(xs([1, end]));
    ylim([min(y,[],"all"),max(y,[],"all")]);
    if nargin == 3, ylim(l), end

    n = size(y, 2);

    for i = 1:20:n
        p.YData = y(:,i);
        drawnow;
    end

end
