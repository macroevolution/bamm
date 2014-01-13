Installation
============

The following instructions assume you are running ``bash`` or a similar
Unix shell. In Linux and Mac OS X systems, ``bash`` is the default shell
(in Mac OS X, you must first open the Terminal application).
For Windows users, we recommend you use `Cygwin <http://www.cygwin.com/>`_,
where ``bash`` is also the default shell.

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
