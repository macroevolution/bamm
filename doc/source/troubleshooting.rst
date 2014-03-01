Troubleshooting
===============
 
Installation
------------

You don't know where to put the bamm directory
..............................................

On OS X, you may move the uncompressed ``bamm`` directory (folder) to your home
directory. For example, if it was uncompressed in your Downloads directory::

    mv ~/Downloads/bamm-1.0.0 ~

To run BAMM within your data directory, you will need to specify its full
location (path)::

    ~/bamm-1.0.0/bamm -c divcontrol.txt

On Windows, you may move the uncompressed "bamm" directory (folder) to your
``C:`` drive. To run BAMM within your data directory, you will need to specify
its full path::

    c:\bamm\bamm -c divcontrol.txt

Running
-------

On OS X, you get an error that says something like "dyld: unknown required load command 0x80000022"
...................................................................................................

This likely means that you are trying to run BAMM on an older version of OS X.
Please run BAMM on OS X 10.5 or greater.

On Windows, you get an error that says "bamm.exe is not a valid Win32 application"
............................................................................................

This likely means that you are trying to run BAMM on an unsupported version
of Windows (such as Windows XP). Please run BAMM on Windows 7 or greater.

You get an error message by BAMM
................................

Most error messages in BAMM are descriptive and therefore it should not be
difficult to correct the problem. If BAMM cannot find the control file or
any file specified within the control file, such as the tree file, make sure
that the file name is spelled correctly and that the location given in the
control file is relative to the directory in which ``bamm`` was called.
Some error messages in BAMM are more cryptic and may require that you contact
us directly.
Please go to `Contact Us <http://bamm-project.org/contact_us.html>`_ on
instructions on how to contact us.

Basic troubleshooting
.....................

The simplest way to localize potential problems is to first try to load your
data without doing anything else. You can do this by setting the following
parameter in your controlfile::

	initializeModel = 0 

This tells BAMM to simply try reading the tree (and associated trait or taxon
sampling datafiles). If this works fine, your problem is probably elsewhere.
The next step is to try initializing your MCMC simulation, without actually
running it. This step involves actually computing likelihoods and priors,
so may reveal problems with those steps::

	initializeModel = 1 
	runMCMC = 0

If these issues fail to resolve the problem, please contact us.
Go to `Contact Us <http://bamm-project.org/contact_us.html>`_ on
instructions on how to contact us.
 
You get an operating system error (e.g., segmentation fault)
............................................................

Other possible errors may arise, not discussed above. This may include errors
that mention "Segmentation fault", or "Floating point error". Please contact
us with information about system errors.
Go to `Contact Us <http://bamm-project.org/contact_us.html>`_ on
instructions on how to contact us.

Analysis
--------

**The MCMC does not converge.** [Under construction.]
