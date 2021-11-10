#pragma once

#include <stdbool.h>
#include <stdint.h>

bool solveStep(
        intmax_t n,
        double const a1[n],
        double const b1[n],
        double a2[n],
        double b2[n],
        double x0,
        double dx,
        double dt);

intmax_t solve(
        intmax_t m,
        intmax_t n,
        double const a0[m],
        double const b0[m],
        double a[n][m],
        double b[n][m],
        double x0,
        double dx,
        double dt);
