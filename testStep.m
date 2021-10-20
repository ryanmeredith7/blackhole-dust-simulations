% Paramaters to set up the initial conditions and solver.
p0 = 0.02;
x0 = 4;
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
a1 = (x0 / a0 ./ xm) .^ 2;
b1 = -asin(sqrt(8*pi*p0/3*x0^3 ./ xm .^ 3 - a1));

a1(xm < x0) = a0 ^ -2;
b1(xm < x0) = -asin(sqrt(8*pi*p0/3 - a0 ^ -2));

plotf = @(a,b) a;

d = plotf(a1, b1);
p = plot(xm, d);
xlim([0, a0]);
ylim([0, 0.2]);

for i = 1:n

    [a2,b2] = solverStep(a1, b1, x0, dx, dt, x, xo2);

    try
        p.YData = plotf(a2, b2);
        drawnow limitrate;
    catch e
        if e.identifier == "MATLAB:class:InvalidHandle"
            warning("Not able to update plot, assumed that user ended early.");
            break
        else
            rethrow(e);
        end
    end

    a1 = a2;
    b1 = b2;

end
