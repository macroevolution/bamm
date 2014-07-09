:orphan:

Download BAMM
=============

Installation instructions are in `Setting Up BAMM <settingup.html>`_.

OS X
----

`bamm-2.0.0-MacOSX.tar.gz
<http://www-personal.umich.edu/~carlosja/bamm-2.0.0-MacOSX.tar.gz>`_

Requires OS X 10.5 or greater.

Windows
-------

`bamm-2.0.0-Windows.zip
<http://www-personal.umich.edu/~carlosja/bamm-2.0.0-Windows.zip>`_

Requires Windows 7 or greater (does not run on Windows XP).

Examples
--------

The following tar.gz and zip files (for compatibility) contain
example data files for diversification and phenotypic evolution analysis.

`bamm-examples.tar.gz
<http://www-personal.umich.edu/~carlosja/bamm-examples.tar.gz>`_

`bamm-examples.zip
<http://www-personal.umich.edu/~carlosja/bamm-examples.zip>`_

Download BAMMtools
==================

BAMMtools will soon be available on CRAN. For now, you can download the full
package here: `BAMMtools_2.0.0.tar.gz
<http://www-personal.umich.edu/~carlosja/BAMMtools_2.0.0.tar.gz>`_
or (for Windows users) `BAMMtools_2.0.0.zip
<http://www-personal.umich.edu/~carlosja/BAMMtools_2.0.0.zip>`_.

You can install this in R in several ways. On Mac OS X, you can open the
Terminal application, navigate to a directory containing the
BAMMtools_2.0.0.tar.gz file, and type "R CMD INSTALL BAMMtools_2.0.0.tar.gz".
Or you can launch the R GUI and go to Packages & Data -> Package
Installer -> Local Source Package, and use the dialog box to navigate to
the BAMMtools_2.0.0.tar.gz file. It should now behave like any other package
you've installed in R, and you can load it with "library(BAMMtools)".

On Windows, first you need to uncompress the zip file
so that you obtain a single directory called "BAMMtools".
Then, download and install
`Rtools <http://cran.r-project.org/bin/windows/Rtools>`_ using RGui.
Next, install and load the R package "devtools" in RGui.
Finally, run ``install(choose.dir())`` and choose the "BAMMtools" directory.
To start using BAMMtools, run ``library(BAMMtools)`` as with any other package.

Download BAMM Source Files
==========================

Access the development source files from our
`GitHub page <https://github.com/macroevolution/bamm>`_.

BAMM Changes
============

Read important changes to BAMM since its release:
`BAMM Changes <changes.html>`_

Previous Versions
=================

You may download previous versions of the binaries and examples
for BAMM and the source code for BAMMtools using the following links:

`bamm-1.0.0-MacOSX.tar.gz
<http://www-personal.umich.edu/~carlosja/bamm-1.0.0-MacOSX.tar.gz>`_

`bamm-1.0.0-Windows.zip
<http://www-personal.umich.edu/~carlosja/bamm-1.0.0-Windows.zip>`_

`bamm-examples-1.0.0.tar.gz
<http://www-personal.umich.edu/~carlosja/bamm-examples-1.0.0.tar.gz>`_

`bamm-examples-1.0.0.zip
<http://www-personal.umich.edu/~carlosja/bamm-examples-1.0.0.zip>`_

`BAMMtools_1.0.1.tar.gz
<http://www-personal.umich.edu/~carlosja/BAMMtools_1.0.1.tar.gz>`_
