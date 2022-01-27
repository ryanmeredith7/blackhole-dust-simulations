function [t1,t2] = timer(a, b, x0, dx, dt)

    arguments

        % Input alpha vector.
        a (:,1) {mustBeNonempty, mustBeReal, mustBeFinite}

        % Input beta vector.
        b (:,1) {mustBeNonempty, mustBeReal, mustBeFinite}

        % Start of the space domain for the solution.
        x0 (1,1) {mustBeNonnegative}

        % Length of the cells in the space direction.
        dx (1,1) {mustBePositive}

        % Time step between steps of the solution.
        dt (1,1) {mustBePositive}

    end

    n = length(a);
    assert(length(b) == n, "Input vectors must be the same length");

    x = (0:dx:(n*dx)).';
    xo2 = x ./ 2;
    xm = x(1:n) + dx/2;

    t1 = NaN;
    t2 = NaN;

    t = 0;

    while isnan(t2)

        theta = 1 - xm .^ 2 .* (a + sin(2 .* b) .^ 2 ./ 4);
        if isnan(t1)
            if any(theta <= 0) && any(theta >= 0)
                t1 = t;
            end
        else
            if all(theta > 0) || all(theta < 0)
                t2 = t;
            end
        end

        [a,b] = solverStep(a, b, x0, dx, dt, x, xo2);
        t = t + dt;

    end

end
