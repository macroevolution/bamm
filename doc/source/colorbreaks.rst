

Color ramps and phylorate plots
===============================


When plotting instantaneous rates along the branches of a phylogeny with ``plot.bammdata`` (as well as with any other BAMMtools function that creates a phylorate plot), it is important to keep in mind how rates are mapped to colors, and what this can mean for interpretation. 

``plot.bammdata`` has a number of options available to bin rates into particular colors (the breaksmethod options). Depending on the dataset, these different options can lead to drastically different phylorate plots. These breaksmethod options are:

1. **linear method**:
Rates across the phylogeny are divided into equal-width bins, and a color from the defined color ramp is associated with each of these bins. The range of this linear color mapping is defined by the fastest and the slowest rates across the tree, regardless of the frequency of these rates. This method can lead to misleading interpretations if the fastest or slowest rates are only represented by a small fraction of the tree, giving the appearance of a lack of rate heterogeneity. 

2. **quantile method**:
Rates are divided into bins, such that each bin has the same density of rates. This method can potentially be misleading in that it might over-emphasize minor rate variation. 


3. **Jenks natural breaks method**:
Rates are divided into bins such that the variance within bins is minimized, while the variance between bins is maximized. 

It is important to understand that while these different color mapping methods can lead to different phylorate plots, the underlying raw data remains the same and any downstream quantitative analyses of rates are completely unaffected by the breaks method chosen. 


We illustrate these differences with an example. This is a BAMM analysis of primate body mass (Vos and Moors 2008). For each breaks method, we generated a histogram showing the frequency of rate values in the dataset, colored according to the selected breaks method. Vertical bars show how the rates have been binned. 

.. _breaksTest:
.. figure:: figs/breaksTest.png
	:width: 700
	:align: center

Now we show how these different methods produce different phylorate figures. 

.. _breaksTestTrees:
.. figure:: figs/breaksTestTrees.png
	:width: 700
	:align: center


How to choose?
..............
The goal is to find the right balance between a visually appealing figure and one that properly portrays the information. We do not want to hide real rate variation, but we also do not want to exaggerate rate heterogeneity "noise". 
One recommendation would be to use the locations of core shifts as a guide for what **should** be emphasized with the colors. In the above figure, the linear method fails to highlight existing rate shifts, and the quantile method exaggerates the existing rate variation. In our tests, the jenks method appears to do a good job in most cases. 






