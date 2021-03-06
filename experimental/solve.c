#include <stdint.h>
#include "mex.h"
#include "libsolve.c"

// MATLAB entry point.
void mexFunction(int nout, mxArray *out[], int nin, const mxArray *in[]) {


    // Checks that the function was called with the correct number of inputs and
    // outputs.
    if (nin != 6) {
        mexErrMsgIdAndTxt("PDESolve:ArgumentError",
                "Requires 6 input arguments");
    }

    if (nout != 2) {
        mexErrMsgIdAndTxt("PDESolve:OutputError", "Requires 2 outputs");
    }

    // Validates and obtains the first two arguments.
    if (
            !mxIsDouble(in[0]) || mxIsEmpty(in[0]) || mxIsComplex(in[0]) ||
            !mxIsDouble(in[1]) || mxIsEmpty(in[1]) || mxIsComplex(in[1])
       ) {
        mexErrMsgIdAndTxt("PDESolve:ArgumentError",
                "The first two arguments must be real non-empty double arrays");
    }

    if (
            mxGetNumberOfDimensions(in[0]) != 2 ||
            mxGetNumberOfDimensions(in[1]) != 2 ||
            mxGetDimensions(in[0])[1] != 1 ||
            mxGetDimensions(in[1])[1] != 1
       ) {
        mexErrMsgIdAndTxt("PDESolve:ArgumentError",
                "The First two arguments must be column vectors");
    }

    uintmax_t m = mxGetNumberOfElements(in[0]);
    if (mxGetNumberOfElements(in[1]) != m) {
        mexErrMsgIdAndTxt("PDESolve:ArgumentError",
                "The initial data vectors must be the same length");
    }

    double *ai = mxGetDoubles(in[0]);
    double *bi = mxGetDoubles(in[1]);

    for (uintmax_t i = 0; i < m; ++i) {
        if (!mxIsFinite(ai[i]) || !mxIsFinite(bi[i])) {
            mexErrMsgIdAndTxt("PDESolve:ArgumentError",
                    "Initial data vectors must be finite");
        }
    }

    // Validates and obtains the next 3 arguments.
    if (
            !mxIsNumeric(in[2]) || !mxIsScalar(in[2]) || mxIsComplex(in[2]) ||
            !mxIsNumeric(in[3]) || !mxIsScalar(in[3]) || mxIsComplex(in[3]) ||
            !mxIsNumeric(in[4]) || !mxIsScalar(in[4]) || mxIsComplex(in[4])
       ) {
        mexErrMsgIdAndTxt("PDESolve:ArgumentError",
                "Arguments 3, 4, and 5 must be real scalar numerical values");
    }

    double x0 = mxGetScalar(in[2]);
    double dx = mxGetScalar(in[3]);
    double dt = mxGetScalar(in[4]);

    if (x0 < 0) {
        mexErrMsgIdAndTxt("PDESolve:ArgumentError",
                "Cannot solve at negative x");
    }

    if (dx <= 0 || dt <= 0) {
        mexErrMsgIdAndTxt("PDESolve:ArgumentError",
                "Steps in time and space must be positive");
    }

    // Validates and obtains the last argument.
    if (!mxIsUint64(in[5]) || !mxIsScalar(in[5]) || mxIsComplex(in[5])) {
        mexErrMsgIdAndTxt("PDESolve:ARgumentError",
                "Argument 6 must be a real scalar uint64, "
                "use uint64(n) instead of just n");
    }

    uintmax_t n = *mxGetUint64s(in[5]);

    if (n < 1) {
        mexErrMsgIdAndTxt("PDESolve:ArgumentError",
                "Number of time steps must be positive");
    }

    // Allocates 2D arrays to hold the results of solve.
    double (*a1)[m] = mxMalloc(n * sizeof(double[m]));
    double (*b1)[m] = mxMalloc(n * sizeof(double[m]));

    // Applies solve, this is the main computation.
    uintmax_t k = solve(m, n, ai, bi, a1, b1, x0, dx, dt);

    // Warns user if the requested number of steps was not completed.
    if (k != n) {
        mexWarnMsgIdAndTxt("PDESolve:ImaginaryNumber",
                "Ran into the imaginary number while computing the boundary "
                "conditions, truncating result");
    }

    // Creates MATLAB arrays to put the output into.
    out[0] = mxCreateDoubleMatrix(m, k, mxREAL);
    out[1] = mxCreateDoubleMatrix(m, k, mxREAL);

    // Gets access to the MATLAB arrays in the form of 1D C arrays.
    double *a2 = mxGetDoubles(out[0]);
    double *b2 = mxGetDoubles(out[1]);

    // Puts the output into the MATLAB arrays.
    for (uintmax_t i = 0; i < k; ++i) {
        for (uintmax_t j = 0; j < m; ++j) {
            a2[i * m + j] = a1[i][j];
            b2[i * m + j] = b1[i][j];
        }
    }

    // Frees the no longer needed 2D arrays.
    mxFree(a1);
    mxFree(b1);

}
