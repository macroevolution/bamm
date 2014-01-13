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
you have CMake installed. Next, in the top directory, run:

    make

This will create a new `build` directory and automatically run `cmake`
to compile BAMM. The executable will be named `bamm`.
You can run `bamm` from this directory, but you may wish to install it
in a more permanent location. To do this, run:

    sudo make install

You may now run `bamm` from any directory in your system.
