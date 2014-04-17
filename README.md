bamm
====

A program for multimodel inference on speciation and trait evolution.
Please see the project's website
([http://bamm-project.org](http://bamm-project.org))
for download options and documentation.

Installation
------------

First, verify that you have [CMake](http://www.cmake.org) installed:

    which cmake

If the location of the CMake program appears (e.g., `/usr/bin/cmake`),
you have CMake installed.

If you would like OpenMP support (highly recommended for running multiple
Markov chains), you'll need to have installed the OpenMP libraries as well as
a compiler supports that it. GCC 4.7 and above supports OpenMP.
For OpenMP support in Clang, see
[http://clang-omp.github.io](http://clang-omp.github.io).

Next, create a `build` directory and go in it:

    mkdir build
    cd build

To build BAMM, run the following:

    cmake ..
    make -j

If you're using Clang for OpenMP support, you will likely need to tell CMake
to use that version of Clang:

    cmake -DCMAKE_CXX_COMPILER=clang++ ..
    make -j

The final executable will be named `bamm`. You can run `bamm` from this
directory, or you may wish to install it in a more permanent location.
To do this, run the following within the `build` directory:

    sudo make install

You may now run `bamm` from any directory in your system.
