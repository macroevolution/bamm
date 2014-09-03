BAMM
====

A program for multimodel inference on speciation and trait evolution.
Please see the project's website
([http://bamm-project.org](http://bamm-project.org))
for the full documentation.

Requirements
------------

In order to compile BAMM,
you need [CMake](http://www.cmake.org) and a C++11 compiler.
You also need a Unix shell (e.g., `bash`) to run the following commands.

Installation
------------

In the project's root directory,
create a new directory called `build` and go into it:

    mkdir build
    cd build

To compile BAMM, run the following commands:

    cmake ..
    make -j

The final executable will be named `bamm`. You may run `bamm` from this
directory, or you may install it in a more permanent location.
To do this, run the following command within the `build` directory:

    sudo make install

You may now run `bamm` from any directory in your system.
