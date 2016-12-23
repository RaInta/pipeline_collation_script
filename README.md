# pipeline_collation_script
Python script to collate continuous gravitational wave search results and output search band and upper limit XML files.

This is a large, multi-purpose script. It's still being refactored, as it violates a number of DRY principles, making it potentially brittle.

This takes search results and compares these to frequency regions known to have serious data quality issues (from a vetoing procedure known as fscan, which identifies elevated non-stationary spectral lines). It checks that each search job has indeed been submitted and completed without error. Then it looks for the template (a fancy matched filter in N-dimensional parameter space) that has the highest detection statistic ('twoF'). It also produces an XML document of loudest twoF values in a designated frequency region, which is subsequently used to execute Monte-Carlo recovery of synthetic signals ('injections') to determine statistical upper limits, in the case of non-detection.

What is cool about this script, and what leads to a huge simplification in computational complexity over the use of blind hash tables, is the use of the bioinformatics library <a href="https://pypi.python.org/pypi/intervaltree">intervaltree</a>. The veto bands XML and search bands XML files are turned into interval trees, then their intersection is computed in <i>O</i>( Nlog(N) ).

A schematic of how this works:

<img src="https://github.com/RaInta/pipeline_collation_script/raw/master/CollateScript_MergingBands_Schematic.png">

An XML file containing the search bands is compared to another XML file containing what we want in terms of upper limits. These are compared to the regions we know to be bad, contained in a veto bands XML. These are, respectively, read in and converted to intervals, containing information on the job number as a data field.

Then, the intersection of the search and the upper limit band is taken. This is then combined with the veto bands. The result is split by <i>any</i> overlap with veto bands and the overlapping segments are removed. In other words, this region is excised (before reading the respective file, indexed by job number). This splitting process is the most labour-intensive one, costing <i>O</i>(N log(N)) operations, where N is the total number of segments.




Thanks, biologists! 
