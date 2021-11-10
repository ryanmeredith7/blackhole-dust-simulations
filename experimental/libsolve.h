#pragma once

#include <stdbool.h>
#include <stdint.h>

bool solveStep(
        uintmax_t n,
        double const a1[n],
        double const b1[n],
        double a2[n],
        double b2[n],
        double x0,
        double dx,
        double dt);

uintmax_t solve(
        uintmax_t m,
        uintmax_t n,
        double const a0[m],
        double const b0[m],
        double a[n][m],
        double b[n][m],
        double x0,
        double dx,
        double dt);
