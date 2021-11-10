#include <stdbool.h>
#include <stdint.h>
#include <math.h>

#include "libsolve.h"

void speed(double x, double bl, double br, double *sl, double *sr);

void speed(double const x, double const bl, double const br, double *sl,
        double *sr) {

    if (floor(bl / M_PI) < floor(br / M_PI)) {
        *sl = x / 2 * pow(sin(bl), 2) / (bl - br);
        *sr = x / 2 * pow(sin(br), 2) / (br - bl);
    } else if (floor(bl / M_PI + 0.5) > floor(br / M_PI + 0.5)) {
        *sl = x / 2 * pow(cos(bl), 2) / (br - bl);
        *sr = x / 2 * pow(cos(br), 2) / (bl - br);
    } else {
        double s;
        if (bl == br) {
            s = x / 2 * sin(2*bl);
        } else {
            s = x / 2 * (pow(sin(br), 2) - pow(sin(bl), 2)) / (br - bl);
        }
        if (s < 0) {
            *sl = s;
            *sr = 0;
        } else {
            *sl = 0;
            *sr = s;
        }
    }

}

bool solveStep(uintmax_t n, double const a1[n], double const b1[n],
        double a2[n], double b2[n], double const x0, double const dx,
        double const dt) {

    double al = a1[0];
    double ar = a1[0];
    double bl = b1[0];
    double br = b1[0];

    double sl, sr;
    speed(x0, bl, br, &sl, &sr);

    double fa = sr * (ar - al);
    double fb = sr * (br - bl);

    for (uintmax_t i = 1; i < n; ++i) {

        al = a1[i-1];
        ar = a1[i];
        bl = b1[i-1];
        br = b1[i];

        speed(x0 + i * dx, bl, br, &sl, &sr);

        a2[i-1] = al - dt * ((fa + sl*(ar - al))/dx + al*sin(2*bl));
        b2[i-1] = bl - dt * ((fb + sl*(br - bl))/dx + 1.5*pow(sin(bl), 2)
                + al/2);

        fa = sr * (ar - al);
        fb = sr * (br - bl);

    }

    al = a1[n-1];
    bl = b1[n-1];

    double sb2 = pow(sin(bl), 2) - dx/(x0 + n*dx)*(3*pow(sin(bl), 2) + al);
    if (-1 <= sb2 && sb2 <= 1) {
        br = -asin(sqrt(sb2));
    } else {
        return false;
    }
    ar = al*(1 + 2*(br - bl)*sin(2*bl)/(3*pow(sin(bl), 2) + al));

    speed(x0 + n * dx, bl, br, &sl, &sr);

    a2[n-1] = al - dt * ((fa + sl*(ar - al))/dx + al*sin(2*bl));
    b2[n-1] = bl - dt * ((fb + sl*(br - bl))/dx + 1.5*pow(sin(bl), 2) + al/2);

    return true;

}

uintmax_t solve(uintmax_t m, uintmax_t n, double const a0[m],
        double const b0[m], double a[n][m], double b[n][m], double const x0,
        double const dx, double const dt) {

    for (uintmax_t i = 0; i < m; ++i) {
        a[0][i] = a0[i];
        b[0][i] = b0[i];
    }

    for (uintmax_t i = 1; i < n; ++i) {
        if (!solveStep(m, a[i-1], b[i-1], a[i], b[i], x0, dx, dt)) {
            return i;
        }
    }

    return n;

}
