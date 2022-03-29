#include <math.h>
#ifdef MATLAB_MEX_FILE
    #define malloc(n) mxMalloc(n)
    #define free(p) mxFree(p)
#else
    #pragma message "Not using MATLAB's heap"
    #include <stdlib.h>
#endif

#include "libsolve.h"

// Type to hold speeds to both the left and right.
typedef struct {
    double l, r;
} speeds;

/* This function calculates the left and right wave speeds given the location,
 * and left and right beta values.
 *
 * Inputs:
 * x is the location in space
 * bl is the left value of beta
 * br is the right value of beta
 *
 * Outputs:
 * Returns a speeds struct s with speeds to the left and right
 */
speeds speed(double x, double bl, double br);

speeds speed(const double x, const double bl, const double br) {

    speeds s;

    // The outer if check for rarefaction waves and calculates the appropriate
    // speeds if there is a rarefaction wave.
    if (floor(bl / M_PI) < floor(br / M_PI)) {
        s.l = x / 2 * pow(sin(bl), 2) / (bl - br);
        s.r = x / 2 * pow(sin(br), 2) / (br - bl);
    } else if (floor(bl / M_PI + 0.5) > floor(br / M_PI + 0.5)) {
        s.l = x / 2 * pow(cos(bl), 2) / (br - bl);
        s.r = x / 2 * pow(cos(br), 2) / (bl - br);
    } else {

        // In this case the wave is only in one direction, so we calculate that
        // in s.
        register double sp;
        if (bl == br) {
            sp = x / 2 * sin(2*bl);
        } else {
            sp = x / 2 * (pow(sin(br), 2) - pow(sin(bl), 2)) / (br - bl);
        }

        // If the wave speed is negative it moves left, if it is positive it
        // moves right.
        if (sp < 0) {
            s.l = sp;
            s.r = 0;
        } else {
            s.l = 0;
            s.r = sp;
        }

    }

    return s;

}

// Simple function to calculate fluctuations from wave speed and left and right
// value of the wave.
double flux(double speed, double left, double right);

double flux(const double s, const double l, const double r) {
    return s * (r - l);
}

// Function that encapsulates the update formula for alpha.
double updateA(double flux, double a, double b, double dx, double dt);

double updateA(const double f, const double a, const double b,
               const double dx, const double dt) {
    return a - dt * (f / dx + a * sin(2 * b));
}

// Function that encapsulates the update formula for beta.
double updateB(double flux, double a, double b, double dx, double dt);

double updateB(const double f, const double a, const double b,
               const double dx, const double dt) {
    return b - dt * (f / dx + 1.5 * pow(sin(b), 2) + a / 2);
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
    register speeds s = speed(x0, bl, br);

    // Variables to hold the fluctuations to the right.
    register double fa = flux(s.r, al, ar);
    register double fb = flux(s.r, bl, br);

    // Here we loop over each cell interface.
    for (uintmax_t i = 1; i < n; ++i) {

        // Update the left and right values of alpha and beta.
        al = a1[i-1];
        ar = a1[i];
        bl = b1[i-1];
        br = b1[i];

        // Calculate the speeds at this interface.
        s = speed(x0 + i * dx, bl, br);

        // Apply the updates to the cell on the left of the interface using
        // fluctuations to the right from the previous interface and the
        // fluctuations to the left from this interface.
        a2[i-1] = updateA(fa + flux(s.l, al, ar), al, bl, dx, dt);
        b2[i-1] = updateB(fb + flux(s.l, bl, br), al, bl, dx, dt);

        // Store the fluctuations to the right for the next step.
        fa = flux(s.r, al, ar);
        fb = flux(s.r, bl, br);

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
    if (0 <= sb2 && sb2 <= 1) {
        br = -asin(sqrt(sb2));
    } else {
        return false;
    }
    ar = al * (1 + 2*(br - bl)*sin(2*bl)/(3*pow(sin(bl), 2) + al));

    // Now we just do the usual calculation of the wave speeds and update
    // formula.
    s = speed(x0 + n * dx, bl, br);

    a2[n-1] = updateA(fa + flux(s.l, al, ar), al, bl, dx, dt);
    b2[n-1] = updateB(fb + flux(s.l, bl, br), al, bl, dx, dt);

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

time_result time(uintmax_t n, const double a[n], const double b[n],
          const double x0, const double dx, const double dt) {

    time_result out = {
        .t1 = NAN,
        .t2 = NAN,
        .success = NOMEM,
    };

    double *a1 = malloc(n * sizeof(double));
    double *b1 = malloc(n * sizeof(double));
    double *a2 = malloc(n * sizeof(double));
    double *b2 = malloc(n * sizeof(double));

    if (!(a1 || b1 || a2 || b2)) {
        free(a1);
        free(b1);
        free(a2);
        free(b2);
        return out;
    }

    for (uintmax_t i = 0; i < n; ++i) {
        a2[i] = a[i];
        b2[i] = b[i];
    }

    register double t = 0;
    register double theta = 0;
    register bool positive = false;
    register bool negative = false;

    for (;;) {

        for (uintmax_t i = 0; i < n; ++i) {

            a1[i] = a2[i];
            b1[i] = b2[i];

            theta = 1 - pow(x0 + (i + 0.5) * dx, 2) * (a1[i] + pow(sin(2 * b1[i]), 2) / 4);

            positive |= theta >= 0;
            negative |= theta <= 0;

        }

        if (isnan(out.t1)) {
            if (positive & negative) {
                out.t1 = t;
            }
        } else {
            if (positive ^ negative) {
                out.t2 = t;
                out.success = SUCCESS;
                break;
            }
        }

        if (!solveStep(n, a1, b1, a2, b2, x0, dx, dt)) {
            out.success = BCERROR;
            break;
        }

        t += dt;

        positive = false;
        negative = false;

    }

    free(a1);
    free(b1);
    free(a2);
    free(b2);

    return out;

}
