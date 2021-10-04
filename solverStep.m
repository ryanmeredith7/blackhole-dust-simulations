% Solves one step for the PDE.
function [a2, b2] = solverStep(a1, b1, x0, dx, dt, x, xo2)

    % This block can be commented out to improve efficiency if you know that the
    % arguments are correct.
    arguments

        % Input alpha vector.
        a1 (:,1) {mustBeNonempty, mustBeReal, mustBeFinite}

        % Input beta vector.
        b1 (:,1) {mustBeNonempty, mustBeReal, mustBeFinite}

        % Start of the space domain for the solution.
        x0 (1,1) {mustBeNonnegative}

        % Length of the cells in the space direction.
        dx (1,1) {mustBePositive}

        % Time step between steps of the solution.
        dt (1,1) {mustBePositive}

        % Optional inputs that can improve efficiency so they don't have to
        % be calculated on each step.

        % Array of locations of cell boundaries.
        x (:,1) {mustBeNonempty, mustBeNonnegative} ...
            = (x0 : dx : x0 + length(a1) * dx).'

        % x over 2.
        xo2 (:,1) {mustBeNonempty, mustBeNonnegative} = 0.5 .* x

    end

    % Gets the  length of the input vectors and checks that they are all the
    % correct lengths. The assertion can be commented out to improve efficiency
    % if you know the arguments are all the same length.
    m = length(a1);
    assert(length(b1) == m && length(x) == m + 1 && length(xo2) == m + 1, ...
        "PDESolve:Step:ArgumentError", ...
        "Both initial data vectors must be the same length.");

    % First we calculate some things we'll be using a lot and store them in
    % arrays. Again not strictly necessary, but done for convenience and
    % efficiency.

    % sin of beta squared.
    sb2 = sin(b1) .^ 2;

    % These are the boundary conditions at the end of the domain, chosen so
    % that the alpha and beta will be fixed when the waves are moving to
    % the left. This should cause density to be 0 at the end point too.
    sb2end = sb2(m) - dx / x(m+1) * (3 * sb2(m) + a1(m));
    bend = -asin(sqrt(sb2end));
    aend = (1 + 2 * (bend - b1(m)) * sin(2 * b1(m)) ...
        / (3 * sb2(m) + a1(m))) * a1(m);

    % Values of beta to the left and right of each cell interface.
    bl = [b1(1); b1];
    br = [b1; bend];

    % Values of sin of beta squared to the left and right of each cell
    % interface.
    sb2l = [sb2(1); sb2];
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
    sr = zeros(m, 1);
    sl = zeros(m, 1);

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
    a2 = a1 - dt .* ((sr .* diff([a1(1); a1]) ...
        + sl .* diff([a1; aend])) ./ dx + sin(2 .* b1) .* a1);
    b2 = b1 - dt .* ((sr .* diff(bl) + sl .* diff(br)) ./ dx ...
        + 1.5 .* sb2 + 0.5 .* a1);

    % As a final check, we make sure the solution is real and if not we
    % throw an error.
    assert(isreal(a2) && isreal(b2), "PDESolve:Step:ImaginaryNumber", ...
        "Ran into the imaginary number.");

end
