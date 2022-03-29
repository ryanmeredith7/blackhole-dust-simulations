# Simulating Black Holes as Dust Clouds with Loop Quantum Gravity Corrections

## Installation

This repository uses MATLAB's [projects](https://www.mathworks.com/help/matlab/projects.html)
feature, so it should be fairly easy to import and use. The simplest way to import this project into
MATLAB click on New>Project>From Git and paste the link to this repository
(`https://github.com/ryanmeredith7/blackhole-dust-simulations`) into the "Repository Path" field.

![Dropdown Menu](photos/FromGit.png)

![Popup Menu](photos/RepoPath.png)

Note that the "Sandbox" field allows you to change where the project will go on your local machine.
Once you have the project in MATLAB, you can update to future changes by simply clicking the "Pull"
button under the project tab.

![Git Pull](photos/Pull.png)

## Usage

The main code is found in the `lib` folder, these are the functions that will do computations. If
you set up the MATLAB project, then these function should be on MATLAB's path, but if you just
manually downloaded the code you may have to do this yourself. The `tests` folder contains example
usages of the functions in `lib` and should also be on MATLAB's path if you're using the project.
Most or all end user scripts, such as the tests, are also accessible through the "Project Shortcuts"
tab.

![Project Shortcuts](photos/Shortcuts.png)

## Experimental C Functions

There are C functions that speed up execution by up to 300%, they require a C compiler that is
recognizable by MATLAB, you can read more about that
[here](https://www.mathworks.com/support/requirements/supported-compilers.html). To use them you
need to run the build script

    >> build

and then the C function should override the MATLAB functions when used. To confirm this, you can use
the commend `which`, for example to see which function is being used when you call `solve`, run the
command

    >> which solve

if the file has a `.m` extension then it is a MATLAB function and if has some sort of `.mex`
extension then it is a C function. To view more information you could also run

    >> which solve -all

where it will also list shadowed functions. To switch back to the MATLAB function, for now you need
to manually delete the created `.mex` files in the `lib` folder.

If you use the C functions, you will likely need to rerun the `build` script after getting the
latest software with pull to make sure you are using the latest version.
