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

ai = 2/3*p0*x0 ./ xs .* sin(4*pi/x0 .* xs);
ai(xs > x0) = 0;

bi = -asin(sqrt(mi - ai));

% Calls the function that solves the equations.
[a, b] = solver(ai, bi, 0, dx, dt, n);

% Calculates the denssity.
bmid = (b(:,1:end-1) + b(:,2:end)) ./ 2;
p = (diff(b, 1, 2) + diff(a, 1, 2) ./ sin(2 .* bmid)) ./ (-4*pi*dt);

% Plays a short movie of the solution values.
figure(Name="Animation of a");
animate(xs, a);
uiwait(msgbox("Press OK to continue.", "Done plotting a"));

figure(Name="Animation of alpha");
animate(xs, sqrt(1 - xs .^ 2 .* a));
uiwait(msgbox("Press OK to continue.", "Done plotting alpha"));

figure(Name="Animation of beta");
animate(xs, b);
uiwait(msgbox("Press OK to continue.", "Done plotting beta"));

figure(Name="Animation of rho");
animate(xs(2:end) - dx/2, diff(xs .^ 3 .* (sin(b) .^ 2 + a)) ./ (8*pi*dx .* (xs(2:end) - dx/2) .^ 2), [0, 0.3]);
uiwait(msgbox("Press OK to continue.", "Done plotting rho"));

function animate(xs, y, l)

    p = plot(xs, y(:,1));
    xlim(xs([1, end]));
    ylim([min(y,[],"all"),max(y,[],"all")]);
    if nargin == 3, ylim(l), end;

    n = size(y, 2);

    for i = 1:5:n
        try
            p.YData = y(:,i);
        catch e
            if e.identifier == "MATLAB:class:InvalidHandle"
                warning("Not able to update plot, assumed that user ended early.");
                break
            else
                rethrow(e);
            end
        end
        drawnow;
    end

end
