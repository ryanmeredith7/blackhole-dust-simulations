mex experimental/solve.c -R2018a
movefile("solve." + mexext, "lib");

mex experimental/time.c -R2018a
movefile("time." + mexext, "lib");
