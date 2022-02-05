#include <stdint.h>
#include "mex.h"
#include "libsolve.c"

// MATLAB entry point.
void mexFunction(int nout, mxArray *out[], int nin, const mxArray *in[]) {


    // Checks that the function was called with the correct number of inputs and
    // outputs.
    if (nin != 5) {
        mexErrMsgIdAndTxt("PDESolve:Time:ArgumentError",
                "Requires 5 input arguments");
    }

    if (nout != 2) {
        mexErrMsgIdAndTxt("PDESolve:Time:OutputError", "Requires 2 outputs");
    }

    // Validates and obtains the first two arguments.
    if (
            !mxIsDouble(in[0]) || mxIsEmpty(in[0]) || mxIsComplex(in[0]) ||
            !mxIsDouble(in[1]) || mxIsEmpty(in[1]) || mxIsComplex(in[1])
       ) {
        mexErrMsgIdAndTxt("PDESolve:Time:ArgumentError",
                "The first two arguments must be real non-empty double arrays");
    }

    if (
            mxGetNumberOfDimensions(in[0]) != 2 ||
            mxGetNumberOfDimensions(in[1]) != 2 ||
            mxGetDimensions(in[0])[1] != 1 ||
            mxGetDimensions(in[1])[1] != 1
       ) {
        mexErrMsgIdAndTxt("PDESolve:Time:ArgumentError",
                "The First two arguments must be column vectors");
    }

    uintmax_t n = mxGetNumberOfElements(in[0]);
    if (mxGetNumberOfElements(in[1]) != n) {
        mexErrMsgIdAndTxt("PDESolve:Time:ArgumentError",
                "The initial data vectors must be the same length");
    }

    double *a = mxGetDoubles(in[0]);
    double *b = mxGetDoubles(in[1]);

    for (uintmax_t i = 0; i < n; ++i) {
        if (!mxIsFinite(a[i]) || !mxIsFinite(b[i])) {
            mexErrMsgIdAndTxt("PDESolve:Time:ArgumentError",
                    "Initial data vectors must be finite");
        }
    }

    // Validates and obtains the next 3 arguments.
    if (
            !mxIsNumeric(in[2]) || !mxIsScalar(in[2]) || mxIsComplex(in[2]) ||
            !mxIsNumeric(in[3]) || !mxIsScalar(in[3]) || mxIsComplex(in[3]) ||
            !mxIsNumeric(in[4]) || !mxIsScalar(in[4]) || mxIsComplex(in[4])
       ) {
        mexErrMsgIdAndTxt("PDESolve:Time:ArgumentError",
                "Arguments 3, 4, and 5 must be real scalar numerical values");
    }

    double x0 = mxGetScalar(in[2]);
    double dx = mxGetScalar(in[3]);
    double dt = mxGetScalar(in[4]);

    if (x0 < 0) {
        mexErrMsgIdAndTxt("PDESolve:Time:ArgumentError",
                "Cannot solve at negative x");
    }

    if (dx <= 0 || dt <= 0) {
        mexErrMsgIdAndTxt("PDESolve:Time:ArgumentError",
                "Steps in time and space must be positive");
    }

    time_result result = time(n, a, b, x0, dx, dt);

    switch (result.success) {
        case NOMEM:
            mexErrMsgIdAndTxt("PDESolve:Time:OutOfMemory",
                    "Ran out of memory");
            break;
        case BCERROR:
            mexErrMsgIdAndTxt("PDESolve:Time:ImaginaryNumber",
                    "Ran into the imaginary number while computing the boundary "
                    "conditions, truncating result");
            break;
        case SUCCESS:
            out[0] = mxCreateDoubleScalar(result.t1);
            out[1] = mxCreateDoubleScalar(result.t2);
            break;
    }

}
