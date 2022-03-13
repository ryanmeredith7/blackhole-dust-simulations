% Paramaters to set up the initial conditions and solver.
x0 = 40;
p0 = 1e-5;
a0 = 50;
dx = 0.1;
dt = 0.02

% Cell boundries
x = (0:dx:a0).';
xo2 = x ./ 2;
% Cell midpoints
xm = x(2:end) - dx/2;

% Creates the initial condition vectors, first the outside case, the the inside
% case.
a = zeros(size(xm));

n = 100;
step = 1e-6;

ps = (p0:step:p0 + step * (n - 1)).';
ts = NaN(1, n);

for i = 1:n

    b = -asin(sqrt(8*pi*ps(i)/3*x0^3 ./ xm .^ 3));
    b(xm < x0) = -asin(sqrt(8*pi*ps(i)/3));

    [t1,t2] = time(a, b, 0, dx, dt);
    disp([t1,t2]);

    if t1 ~= 0
        ts(i) = t2 - t1;
    end

end

i = isnan(ts);

ps(i) = [];
ts(i) = [];

ms = 4/3*pi*x0^3 .* ps;

f = @(p, x) p(1) .* x .^ 2 + p(2) .* x + p(3);

params = lsqcurvefit(f, [8/3*pi, 0, 0], ms, ts);

disp(params(1))
disp(8*pi/3)
