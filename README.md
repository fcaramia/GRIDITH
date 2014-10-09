Loi Lab @ Peter MacCallum Cancer Centre

GRIDITH
===

Deciphering clonality in aneuploid tumors

Developed by Ingrid Lonnstedt, Franco Caramia 

Publication:
Deciphering clonality in aneuploid tumors using SNP array and sequencing data
http://www.ncbi.nlm.nih.gov/pubmed/25270265



=== Input Data ===

a data.frame consisting of 7 columns:

Chromosome: integer - Chromosome where segment is located

Start.bp : interger - starting position of segment

End.bp : integer - end position of segment

a1: real - minor copy number estimate on arbitrary scale

a2: real - major copy number estimate on arbitrary scale

W: real - length of segment divided by the total genome length

length: integer - length of segment 

=== Functions === 

start.alphamax.f: Enables manual input of starting values for the rotation function optimizeRotation(), which in turn normalizes SNP Array segment CNs for BAF bias (seen as skewness in grid plots).

plot.transformed: Graphical representation of roration starting values from start.alphamax.f().

setWeights: Function which visualizes potential parameter settings for optimizeRotation(), which in turn normalizes SNP Array segment CNs for BAF bias (seen as skewness in grid plots).

optimizeRotation: Function to normalize SNP Array segment CNs for BAF bias (seen as skewness in grid plots). The function requires starting values e.g. as generated via the function start.alphamax.f(). It will optimize the fit of a rectangular lattice pattern to a grid plot by minimizing the summed distance from each grid point to its closest lattice point through Tukey's robust function for M-estimation. See paper... Note that the function can only find a locally optimized rotation. The starting values are important to find the globally best fit.

setCutsUnscaled: Function to set the upper, lower and right side limit of the grid plot rectangle within which a regular grid pattern is evident. The grid plot refers to a plot of minor versus major segment CN plot.

subclonalCNDist: Define a Mahalanobis distance from each lattice point which defines the cutpoint between noise deviation from the lattice point and deviation due to heterozygsity. Choose a cut so were QQ-plot points are clearly no longer on the y=x diagonal.

setType: Function to split segments into either of three types: Segments for which the CNs follow a regular grid pattern (which suggests these segments have CN alteration in one single subclone: the main subclone), segments which break the regular grid pattern (which suggests these segments have CN alteration in not only the main subclone), and segments which are not possible to classify. The function also calculates ITH, and endpoint which estimates the fraction of the genome which has CN alteration in more than one subclone.

