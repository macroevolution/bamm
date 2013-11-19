Running BAMM
============

1.  Speciation-extinction analyses require a plain text control file and a Newick-formatted tree. Trait analyses additionally require a file containing tip values associated with the tree.
   
2. To run speciation-extinction analyses, run the following command::

       ./bamm speciationextinction -control controlfile
       
   Or for phenotypic analyses::

       ./bamm trait -control controlfile

