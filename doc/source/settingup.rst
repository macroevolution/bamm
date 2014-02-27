.. _bammsetup:

Setting Up BAMM
===============

The following instructions assume you are running ``bash`` or a similar
Unix shell. In Linux and Mac OS X systems, ``bash`` is the default shell
(in Mac OS X, you must first open the Terminal application).
For Windows users, we recommend you use `Cygwin <http://www.cygwin.com/>`_,
where ``bash`` is also the default shell.

Installation From Binary
------------------------

Mac OS X
........

1. `Download <http://bamm-project.org/download.html>`_ the Mac OS X
   tar.gz file containing a binary file of BAMM.

2. This file is compressed as a tar.gz file. To uncompress it,
   run the following command from the directory in which BAMM was downloaded
   (the actual version of BAMM may be different)::

       tar -xzf bamm-1.0.0-MacOSX.tar.gz

   This will uncompress the single file ``bamm`` into the current directory.

3. You may copy the file ``bamm`` to any directory in your system.
   See below for instructions on how to run BAMM.

Windows
.......

Binary file for Windows will be available soon.
   
Installation From Source
------------------------

Use this option to compile and install BAMM from its source files.

1. `Download <http://bamm-project.org/download.html>`_ the latest available
   version of BAMM.

2. This file is compressed as a tar.gz file. To uncompress it,
   run the following command from the directory in which BAMM was downloaded
   (the actual version of BAMM may be different)::

       tar -xzf bamm-1.0.0.tar.gz

   This will create the directory bamm-1.0.0.
   
3. Verify that you have `CMake <http://www.cmake.org>`_ installed::

       which cmake

   If the location of the CMake program appears (e.g., ``/usr/bin/cmake``),
   you have CMake installed. You will also need a C++ compiler,
   such as `GCC <http://gcc.gnu.org/>`_ or `LLVM <http://llvm.org/>`_.

4. Go into the top level directory and build the program by running ``make``::

       cd bamm-1.0.0
       make

   This will create the directory ``build`` and automatically run ``cmake``
   to compile BAMM. The executable will be named ``bamm``.

5. You can run ``bamm`` from this directory, or you may wish to install it
   in a more permanent location. To do this, run::

       sudo make install

   You may now run ``bamm`` from any directory in your system. See the
   `Quick-start guide to BAMM <http://bamm-project.org/quickstart.html>`_
   to learn how to configure and run BAMM.

Running
-------

To run ``bamm``, you must always specify a *contol* file. For example,
if your control file is named ``divcontrol.txt``, run the following::

    bamm -c divcontrol.txt

Note that if ``bamm`` is not installed in a common location, you may need
to run ``bamm`` from the directory in which it exists as follows::

    ./bamm -c divcontrol.txt

Any option in the control file may be overridden in the command-line
by prefixing the option name by ``--``, followed by the new value.
For example, to set the seed to 1234::

    ./bamm -c divcontrol.txt --seed 1234

To set the initial lambda at the root of the tree to 0.05
and the print frequency to 5000, use::

    ./bamm -c divcontrol.txt --lambdaInit0 0.05 --printFreq 5000

When run, BAMM produces a file named ``run_info.txt`` that logs
the command-line call used, the random seed, the start and end
time-stamps, and a list of parameters/options and their values.
