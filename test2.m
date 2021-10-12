% Paramaters to set up the initial conditions and solver.
p0 = 0.08;
x0 = 4;
a0 = 5;
dx = 0.01;
dt = 0.001;
n = 5000;

% Midpoints of the cells.
xs = (dx/2:dx:a0).';

% Here we set up the initial conditions.
% these have nothing to do with mass, I just needed to name it something
mi = 8*pi*p0/3*x0^3 ./ xs .^ 3;
mi(xs < x0) = 8*pi*p0/3;

ai = sin(4*pi/a0 .* xs) ./ (a0 .* xs);

bi = -asin(sqrt(mi - ai));

% Calls the function that solves the equations.
[a, b] = solver(ai, bi, 0, dx, dt, n);

% Calculates the denssity.
bmid = (b(:,1:end-1) + b(:,2:end)) ./ 2;
p = (diff(b, 1, 2) + diff(a, 1, 2) ./ sin(2 .* bmid)) ./ (-4*pi*dt);

% Plays a short movie of the solution values.
disp("Plotting a");
animate(xs, a);
disp("Done a");
uiwait(gcf);

disp("Plotting alpha");
animate(xs, sqrt(1 - xs .^ 2 .* a));
disp("Done alpha");
uiwait(gcf);

disp("Plotting beta");
animate(xs, b);
disp("Done beta");
uiwait(gcf);

disp("Plotting rho");
animate(xs(2:end) - dx/2, diff(xs .^ 3 .* (sin(b) .^ 2 + a)) ./ (8*pi*dx .* (xs(2:end) - dx/2) .^ 2), [0, 0.3]);
disp("Done rho");
uiwait(gcf);

disp("Plotting p");
animate(xs, p, [0, 0.3]);
disp("Done p");
uiwait(gcf);

function animate(xs, y, l)

    p = plot(xs, y(:,1));
    xlim(xs([1, end]));
    ylim([min(y,[],"all"),max(y,[],"all")]);
    if nargin == 3, ylim(l), end;

    n = size(y, 2);

    for i = 1:5:n
        p.YData = y(:,i);
        drawnow;
    end

end
