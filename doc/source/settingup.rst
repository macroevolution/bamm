.. _bammsetup:

Setting Up BAMM
===============

The following instructions assume you are running ``bash`` or a similar
Unix shell. In Linux and Mac OS X systems, ``bash`` is the default shell
when you open the Terminal application.
For Windows users, we recommend you use `Cygwin <http://www.cygwin.com/>`_,
where ``bash`` is also the default shell.

Installation From Binary
------------------------

The compressed file you downloaded containing the binary file
contains only the executable program.
We strongly recommend that you also download the examples
to use as templates for the control files necessary to run BAMM.
Go to the `Download <http://bamm-project.org/download.html>`_ page
to download the examples.

OS X
....

1. `Download <http://bamm-project.org/download.html>`_ the Mac OS X
   tar.gz file containing the binary file for BAMM.

2. This file is a compressed tar.gz file. To uncompress it,
   run the following command from the directory in which BAMM was downloaded
   (the actual version of BAMM may be different)::

       tar -xzf bamm-1.0.0-MacOSX.tar.gz

   This will uncompress the single file ``bamm`` into the current directory.

3. You may copy the file ``bamm`` to any directory in your system.
   See below and the `Quick-start guide to BAMM
   <http://bamm-project.org/quickstart.html>`_
   to learn how to configure and run BAMM.

Windows
.......

1. `Download <http://bamm-project.org/download.html>`_ the Windows .zip file
   containing the binary file and required library files for BAMM.

2. This file is a compressed zip file. Extract the contents of the zip file
   and move it to a location you can access from the command-line.

3. Open the command-line in Windows by clicking on the Start button
   and typing "cmd" in the search text box. Use the ``cd`` command
   to go to the directory where BAMM is located.
   See below and the `Quick-start guide to BAMM
   <http://bamm-project.org/quickstart.html>`_
   to learn how to configure and run BAMM.
   
Installation From Source
------------------------

Use this option to compile and install BAMM from its source files.

1. `Download <http://bamm-project.org/download.html>`_ the latest available
   version of the BAMM source files.

2. This file is compressed as a tar.gz file. To uncompress it,
   run the following command from the directory in which BAMM was downloaded
   (the actual version of BAMM may be different)::

       tar -xzf bamm-1.0.0.tar.gz

   This will create the directory bamm-1.0.0.
   
3. Verify that you have `CMake <http://www.cmake.org>`_ installed::

       which cmake

   If the location of the CMake program appears (e.g., ``/usr/bin/cmake``),
   you have CMake installed. You will also need a C++ compiler,
   such as `GCC <http://gcc.gnu.org/>`_ or `Clang <http://clang.llvm.org/>`_.

4. Go into the bamm-1.0.0 directory and create a build directory where
   the final executable will be created, then go into that directory::
   
       mkdir build
       cd build

5. Run ``cmake``, and if successful, run ``make``::

       cmake ..
       make -j

   The ``-j`` option will compile in parallel.

6. You can run ``bamm`` from this directory, or you may wish to install it
   in a more permanent location::

       sudo make install

   See below and the `Quick-start guide to BAMM
   <http://bamm-project.org/quickstart.html>`_
   to learn how to configure and run BAMM.

Running
-------

To run ``bamm``, you must always specify a *contol* file. For example,
if your control file is named ``divcontrol.txt``, run the following::

    bamm -c divcontrol.txt

Note that if ``bamm`` is not installed in a common location, you may need
to run ``bamm`` from the directory in which it exists as follows::

    ./bamm -c divcontrol.txt

Any file names specified in the control file are relative to the directory
in which ``bamm`` was called, which may not be the same location as where
the executable ``bamm`` nor the control file reside.

Any option in the control file may be overridden in the command-line
by prefixing the option name by ``--``, followed by the new value.
For example, to set the seed to 1234, run::

    ./bamm -c divcontrol.txt --seed 1234

To set the initial lambda at the root of the tree to 0.05
and the print frequency to 5000, run::

    ./bamm -c divcontrol.txt --lambdaInit0 0.05 --printFreq 5000

When run, BAMM produces a file named ``run_info.txt`` that logs
the command-line call used, the random seed, the start and end
time-stamps, and a list of parameters/options and their values.
