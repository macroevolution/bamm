.. _bammsetup:

Setting Up BAMM
===============

Check your BAMM version
-----------------------

If you have previously downloaded or installed BAMM,
you may check which version of BAMM you have
using the ``--version`` option in ``bamm``::

    bamm --version

If ``bamm`` is not automatically found by your operating system,
you need to specify the directory in which ``bamm`` is located::

    ~/bamm/build/bamm --version

The previous command assumes that ``bamm`` is located
in the ``~/bamm/build`` directory.
If ``bamm`` is in your current directory, you may run::

    ./bamm --version

Installation from Homebrew (OS X only)
------------------------------------------

#. Install `Homebrew <http://brew.sh>`_, if not installed.

#. Run the following commands::

       brew tap macroevolution/bamm
       brew install bamm

   If you've installed BAMM before
   but you want to upgrade to a new version, run::

       brew update
       brew upgrade bamm

#. You may now run ``bamm`` from any directory in your system.
   See :ref:`running` below and the `Quick-start guide to BAMM
   <http://bamm-project.org/quickstart.html>`_
   to learn how to configure and run BAMM.

Installation from Executable
----------------------------

The compressed file you downloaded contains the executable program only.
We recommend that you also download the example files
to use as templates for the control files necessary to run BAMM.
Go to the `Download <http://bamm-project.org/download.html>`_ page
to download the examples.

OS X
....

#. `Download <http://bamm-project.org/download.html>`_ the Mac OS X
   tar.gz file containing the executable file for BAMM.

#. This file is a compressed tar.gz file. To uncompress it,
   open the Terminal application and run the following command
   from the directory in which BAMM was downloaded
   (the actual version of BAMM may be different)::

       tar -xzf bamm-1.0.0-MacOSX.tar.gz

   This will uncompress the single file ``bamm`` into the current directory.

#. You may copy or move the file ``bamm`` to any directory in your system.
   See :ref:`running` below and the `Quick-start guide to BAMM
   <http://bamm-project.org/quickstart.html>`_
   to learn how to configure and run BAMM.

**Note:** Test that ``bamm`` will run on your system by simply running
``./bamm`` in the directory in which the program is located.
If it prints out a message on "Usage", the BAMM program is working correctly.

Windows
.......

#. `Download <http://bamm-project.org/download.html>`_ the Windows .zip file
   containing the executable file and required library files for BAMM.

#. This file is a compressed zip file. Extract the contents of the zip file
   and move them to a location you can access from the command-line.
   Note that the DLL files must present together with the ``bamm.exe`` file.

#. Open the command-line program in Windows by clicking on the Start button
   and typing "cmd" in the search text box. Use the ``cd`` command
   to go to the directory where BAMM is located.
   See :ref:`running` below and the `Quick-start guide to BAMM
   <http://bamm-project.org/quickstart.html>`_
   to learn how to configure and run BAMM.

**Note:** Test that ``bamm`` will run on your system by simply running
``bamm`` in the directory in which the program is located.
If it prints out a message on "Usage", the BAMM program is working correctly.

Installation from Source
------------------------

The following instructions assume you are running ``bash``
or a similar Unix shell.
In Linux and Mac OS X systems, ``bash`` is the default shell
when you open the Terminal application.
We also assume you have `Git <http://git-scm.com>`_ and
`CMake <http://www.cmake.org>`_ installed.

#. Get the latest version of BAMM::

       git clone https://github.com/macroevolution/bamm.git

   This will create the directory ``bamm`` in your current directory.
   
#. Go into the ``bamm`` directory and create a build directory where
   the final executable will be created, then go into that directory::
   
       mkdir build
       cd build

#. Run ``cmake``, and if successful, run ``make``::

       cmake ..
       make -j

   The ``-j`` option will compile in parallel.

#. You can run ``bamm`` from this directory, or you may wish to install it
   in a more permanent location::

       sudo make install

   You may now run ``bamm`` from any directory in your system.
   See :ref:`running` below and the `Quick-start guide to BAMM
   <http://bamm-project.org/quickstart.html>`_
   to learn how to configure and run BAMM.

.. _running:

Running
-------

To run ``bamm``, you should specify a *contol* file. For example,
if your control file is named ``divcontrol.txt``, run the following::

    bamm -c divcontrol.txt

Note that if ``bamm`` is not installed in a common location,
you need to specify the directory in which ``bamm`` is located::

    ~/bamm/build/bamm -c divcontrol.txt

The previous command assumes that ``bamm`` is located
in the ``~/bamm/build`` directory.
If ``bamm`` is in your current directory, you may run::

    ./bamm -c divcontrol.txt

Any file names specified in the control file are relative to the directory
in which ``bamm`` was called, which may not be the same location where
the executable ``bamm`` or the control file are located.

Any option in the control file may be overridden in the command-line
by prefixing the option name by ``--``, followed by the new value.
For example, to set the seed to 1234, run::

    bamm -c divcontrol.txt --seed 1234

To set the initial lambda at the root of the tree to 0.05
and the print frequency to 5000, run::

    bamm -c divcontrol.txt --lambdaInit0 0.05 --printFreq 5000

When run, BAMM produces a file named ``run_info.txt`` that logs
the command-line call used, the random seed, the start and end
time-stamps, and a list of parameters/options and their values.
