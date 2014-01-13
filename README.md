bamm
====

A program for multimodel inference on speciation and trait evolution.
Please see the project's website
([http://www.bamm-project.org](http://www.bamm-project.org))
for download options and documentation.

Installation
------------

First, verify that you have [CMake](http://www.cmake.org) installed:

    which cmake

If the location of the cmake program appears, you have CMake installed.
Next, in the top directory, run:

    make

This will create a new `build` directory and run `cmake` to compile `bamm`.
You can run bamm from this directory, but you may wish to install it
on a more permanent location. To do this, run:

    sudo make install

You may now run bamm from any directory in your system.
