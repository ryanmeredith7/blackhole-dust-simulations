% This function solves a PDE relating to black holes.
function [a,b] = solver(ai, bi, x0, dx, dt, n)

    arguments

        % Initial alpha vector.
        ai (:,1) {mustBeNonempty, mustBeReal, mustBeFinite}

        % Initial beta vector.
        bi (:,1) {mustBeNonempty, mustBeReal, mustBeFinite}

        % Start of the space domain for the solution.
        x0 (1,1) {mustBeNonnegative}

        % Length of the cells in the space direction.
        dx (1,1) {mustBePositive}

        % Time step between steps of the solution.
        dt (1,1) {mustBePositive}

        % Desired number of time steps to preform.
        n (1,1) {mustBeInteger, mustBePositive}

    end

    % Gets the length of the initial vectors and makes sure they are the same.
    m = length(ai);
    assert(length(bi) == m, "PDESolve:ArgumentError", ...
        "Both initial data vectors must be the same length.");

    % Creates the arrays for the solution.
    a = zeros(m, n);
    b = zeros(m, n);
    a(:,1) = ai;
    b(:,1) = bi;

    % Array of locations of the cell interfaces, and those numbers over 2. Not
    % necessary to calculate "x over 2" here, done for efficiency reasons.
    x = (x0 : dx : x0 + m * dx).';
    xo2 = 0.5 .* x;

    for i = 1:n-1

        % First we calculate some things we'll be using a lot and store them in
        % arrays. Again not strictly necessary, but done for convenience and
        % efficiency.

        % sin of beta squared.
        sb2 = sin(b(:,i)) .^ 2;

        % These are the boundary conditions at the end of the domain, chosen so
        % that the alpha and beta will be fixed when the waves are moving to
        % the left. This should cause density to be 0 at the end point too.
        sb2end = sb2(m) - dx / x(m+1) * (3 * sb2(m) + a(m,i));
        bend = -asin(sqrt(sb2end));
        aend = (1 + 2 * (bend - b(m,i)) * sin(2 * b(m,i)) ...
            / (3 * sb2(m) + a(m,i))) * a(m,i);

        % Boundary conditions on the left, set so that a and b are constant in
        % space.
        a0 = a(1,i);
        b0 = b(1,i);

        % Values of beta to the left and right of each cell interface.
        bl = [b0; b(:,i)];
        br = [b(:,i); bend];

        % Values of sin of beta squared to the left and right of each cell
        % interface.
        sb2l = [sin(b0)^2; sb2];
        sb2r = [sb2; sb2end];

        % Next, we calculate the wave speeds.

        % Most will be given by this expression.
        s = xo2 .* (sb2r - sb2l) ./ (br - bl);

        % Fix the speeds where bl = br as they would be NaN otherwise.
        j = bl == br;
        s(j) = xo2(j) .* sin(2 .* bl(j));

        % Now we calculate the speeds of the left going waves and right going
        % waves.

        % First we split the previous waves speeds.
        sr = zeros(m + 1, 1);
        sl = zeros(m + 1, 1);

        % Positive wave speeds are moving to the right.
        j = s > 0;
        sr(j) = s(j);

        % Negative wave speeds are moving to the left.
        j = s < 0;
        sl(j) = s(j);

        % Next we apply fixes for transonic rarefaction waves, which have
        % components moving to both the left and right. It can be shown that
        % this only happens when one of two conditions are true.

        % This condition is necessary and sufficient for bl to be less than br
        % and the minimum of sin squared beta between bl and br to be 0, so we
        % split the speeds with a mid point of 0.
        j = floor(bl ./ pi) < floor(br ./ pi);
        sr(j) = xo2(j) .* sb2r(j) ./ (br(j) - bl(j));
        sl(j) = xo2(j) .* sb2l(j) ./ (bl(j) - br(j));

        % This condition is necessary and sufficient for bl to be greater than
        % br and the maximum of sin squared beta between bl and br to be 1, so
        % we split the speeds with a mid point of 1.
        j = floor(bl ./ pi + 0.5) > floor(br ./ pi + 0.5);
        sr(j) = xo2(j) .* (1 - sb2r(j)) ./ (bl(j) - br(j));
        sl(j) = xo2(j) .* (1 - sb2l(j)) ./ (br(j) - bl(j));

        % Next we remove waves moving off the end of the domain.
        sr(m+1) = [];
        sl(1) = [];

        % Finally, we simply use these speeds in our update formula.
        a(:,i+1) = a(:,i) - dt .* ((sr .* diff([a0; a(:,i)]) ...
            + sl .* diff([a(:,i); aend])) ./ dx ...
            + sin(2 .* b(:,i)) .* a(:,i));
        b(:,i+1) = b(:,i) - dt .* ((sr .* diff(bl) + sl .* diff(br)) ./ dx ...
            + 1.5 .* sb2 + 0.5 .* a(:,i));

        % As a final check, we make sure the solution is real and if not we
        % truncate the output array and end early, warning the user.
        if ~(isreal(b(:,i+1)) && isreal(a(:,i+1)))
            a = a(:,1:i);
            b = b(:,1:i);
            warning("PDESolve:ImaginaryNumber", ...
                "Ran into the imaginary number, ending early.");
            break
        end

    end

end
