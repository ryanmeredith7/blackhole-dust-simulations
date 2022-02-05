#pragma once

#include <stdbool.h>
#include <stdint.h>

typedef struct {
    enum {SUCCESS, NOMEM, BCERROR} success;
    double t1, t2;
} time_result;

/* This function solves a step of a solution to the PDE.
 *
 * Inputs:
 *
 * n is the length of the input and output arrays, all of which must be the same
 * length
 *
 * a1 is the initial vector of alpha, of length n
 *
 * b1 is the initial vector of beta, of length n
 *
 * x0 is the left side of the solution range in space, must be non-negative
 *
 * dx is the step in space, that is the cell distance, must be positive
 *
 * dt is the step in time, must be positive
 *
 * Outputs:
 *
 * a2 must be an array of length n that the updated alpha will be put into
 *
 * b2 must be an array of length n that the updated beta will be put into
 *
 * Returns true if the function was successful and false if it wasn't.
 */
bool solveStep(
        uintmax_t n,
        const double a1[n],
        const double b1[n],
        double a2[n],
        double b2[n],
        double x0,
        double dx,
        double dt);

/* This function will solve the PDE for n steps.
 *
 * Inputs:
 *
 * m is the length of the input arrays and the inner dimension of the 2D output
 * arrays, which must all have the same length
 *
 * n in the number of time steps to compute, which must also be the outer
 * dimension of the output arrays
 *
 * ai is the initial vector of alpha, of length n
 *
 * bi is the initial vector of beta, of length n
 *
 * x0 is the left side of the solution range in space, must be non-negative
 *
 * dx is the step in space, that is the cell distance, must be positive
 *
 * dt is the step in time, must be positive
 *
 * Outputs:
 *
 * a must be an m by n 2D array that alpha will be put into
 *
 * b must be an m by n 2D array that beta will be put into
 *
 * Returns the number of suuccessfuly completed steps.
 */
uintmax_t solve(
        uintmax_t m,
        uintmax_t n,
        const double a0[m],
        const double b0[m],
        double a[n][m],
        double b[n][m],
        double x0,
        double dx,
        double dt);

time_result time(
        uintmax_t n,
        const double a[n],
        const double b[n],
        double x0,
        double dx,
        double dt);
