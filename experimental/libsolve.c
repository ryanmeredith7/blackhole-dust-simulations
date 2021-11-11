#include <stdbool.h>
#include <stdint.h>
#include <math.h>

#include "libsolve.h"

/* This function calculates the left and right wave speeds given the location,
 * and left and right beta values.
 *
 * Inputs:
 * x is the location in space
 * bl is the left value of beta
 * br is the right value of beta
 *
 * Outputs:
 * sl will contain the speed to the left
 * sr will contain the speed to the right
 */
void speed(const double x, const double bl, const double br,
           double *restrict sl, double *restrict sr);

void speed(const double x, const double bl, const double br,
           double *restrict sl, double *restrict sr) {

    // The outer if check for rarefaction waves and calculates the appropriate
    // speeds if there is a rarefaction wave.
    if (floor(bl / M_PI) < floor(br / M_PI)) {
        *sl = x / 2 * pow(sin(bl), 2) / (bl - br);
        *sr = x / 2 * pow(sin(br), 2) / (br - bl);
    } else if (floor(bl / M_PI + 0.5) > floor(br / M_PI + 0.5)) {
        *sl = x / 2 * pow(cos(bl), 2) / (br - bl);
        *sr = x / 2 * pow(cos(br), 2) / (bl - br);
    } else {

        // In this case the wave is only in one direction, so we calculate that
        // in s.
        register double s;
        if (bl == br) {
            s = x / 2 * sin(2*bl);
        } else {
            s = x / 2 * (pow(sin(br), 2) - pow(sin(bl), 2)) / (br - bl);
        }

        // If the wave speed is negative it moves left, if it is positive it
        // moves right.
        if (s < 0) {
            *sl = s;
            *sr = 0;
        } else {
            *sl = 0;
            *sr = s;
        }

    }

}

// One step of solution, see interface file for more details.
bool solveStep(uintmax_t n, const double a1[n], const double b1[n],
               double a2[n], double b2[n], const double x0, const double dx,
               const double dt) {

    // Variables that hold the left and right values of alpha and beta. At the
    // left boundary we make both constant in space.
    register double al = a1[0];
    register double ar = a1[0];
    register double bl = b1[0];
    register double br = b1[0];

    // Variables to hold the wave speeds to the left and right.
    double sl, sr;
    speed(x0, bl, br, &sl, &sr);

    // Variables to hold the fluctuations to the right.
    register double fa = sr * (ar - al);
    register double fb = sr * (br - bl);

    // Here we loop over each cell interface.
    for (uintmax_t i = 1; i < n; ++i) {

        // Update the left and right values of alpha and beta.
        al = a1[i-1];
        ar = a1[i];
        bl = b1[i-1];
        br = b1[i];

        // Calculate the speeds at this interface.
        speed(x0 + i * dx, bl, br, &sl, &sr);

        // Apply the updates to the cell on the left of the interface using
        // fluctuations to the right from the previous interface and the
        // fluctuations to the left from this interface.
        a2[i-1] = al - dt * ((fa + sl*(ar - al))/dx + al*sin(2*bl));
        b2[i-1] = bl - dt * ((fb + sl*(br - bl))/dx + 1.5*pow(sin(bl), 2)
                             + al/2);

        // Store the fluctuations to the right for the next step.
        fa = sr * (ar - al);
        fb = sr * (br - bl);

    }

    // The last interface requires more care because of the boundary.
    al = a1[n-1];
    bl = b1[n-1];

    // Here we compute the right boundary conditions. There were chosen so that
    // alpha and beta will be constant at the end points when the waves are
    // moving to the left. sb2 is sin(b)^2 at the boundary, we need to check
    // that this value makes sense before continuing.
    register double sb2 = pow(sin(bl), 2)
                          - dx / (x0 + n*dx) * (3*pow(sin(bl), 2) + al);
    if (-1 <= sb2 && sb2 <= 1) {
        br = -asin(sqrt(sb2));
    } else {
        return false;
    }
    ar = al * (1 + 2*(br - bl)*sin(2*bl)/(3*pow(sin(bl), 2) + al));

    // Now we just do the usual calculation of the wave speeds and update
    // formula.
    speed(x0 + n * dx, bl, br, &sl, &sr);

    a2[n-1] = al - dt * ((fa + sl*(ar - al))/dx + al*sin(2*bl));
    b2[n-1] = bl - dt * ((fb + sl*(br - bl))/dx + 1.5*pow(sin(bl), 2) + al/2);

    return true;

}

// Applies solveStep n times, see interface file for more details.
uintmax_t solve(uintmax_t m, uintmax_t n, const double a0[m],
                const double b0[m], double a[n][m], double b[n][m],
                const double x0, const double dx, const double dt) {

    // Set the first array in the output to the initial data.
    for (uintmax_t i = 0; i < m; ++i) {
        a[0][i] = a0[i];
        b[0][i] = b0[i];
    }

    // Repeatedly applies solveStep, n times, checking to make sure it was
    // successful. Returns the number of valid (successful) steps completed.
    for (uintmax_t i = 1; i < n; ++i) {
        if (!solveStep(m, a[i-1], b[i-1], a[i], b[i], x0, dx, dt)) {
            return i;
        }
    }

    return n;

}
