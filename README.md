BAMM
====

A program for multimodel inference on speciation and trait evolution.
Please see the project's website
([http://bamm-project.org](http://bamm-project.org))
for download options and documentation.

Requirements
------------

In order to compile BAMM, you need [CMake](http://www.cmake.org)
and a C++11 compiler.

Installation
------------

Create a `build` directory and go in it:

    mkdir build
    cd build

To build BAMM, run the following:

    cmake ..
    make -j

The final executable will be named `bamm`. You can run `bamm` from this
directory, or you may wish to install it in a more permanent location.
To do this, run the following within the `build` directory:

    sudo make install

You may now run `bamm` from any directory in your system.
