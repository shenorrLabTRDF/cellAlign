# cellAlign
global and local alignment of single cell trajectories

author: "[Lindsay Moore, Ayelet Alpert](http://shenorrlab.technion.ac.il/)"

date: "January 12, 2018"

---
To use this package, you 
will need the R statistical computing environment (version 3.0 or later)
and several packages available through Bioconductor and CRAN.

This release supports Windows 10. cellAlign has not been tested with Mac and Linux operating systems.

---

# Abstract 
High-dimensional measurements of single cells such as RNA sequencing or mass cytometry are often used to generate detailed trajectories of complex biological processes from heterogeneous cell populations. The high-dimensional information enables ordering of single cells based on continuous changes in measured mRNA or protein expression.  CellAlign is a tool for quantitative comparison of expression dynamics within or between single-cell trajectories. The **input** to the CellAlign workflow is any **trajectory vector** that orders single cell expression with a pseudo-time spacing and the **expression matrix** for the cells used to define the trajectory. This vignette provides an overview of a single cell RNA-Seq analysis workflow with CellAlign describing both whole trajectory and partial trajectory alignments. 

# To begin
cellAlign 
First you will need to download the package from github 

Start by installing the devtools package from CRAN and load it

Installation should take less than one minute.

```{r eval=FALSE}
install.packages("devtools")
library(devtools)
```

install cellAlign and load it

```{r eval=FALSE}
install_github("lindsaysmoore/cellAlign")
```
```{r}
library(cellAlign)

```

# Introduction 
Cell Align is a package that takes as input a scaled trajectory vector that orders a gene expression matrix along an arbitrary pseudotime scaled to [0,1], and the gene expression matrix. cellAlign has 3 essential steps:

1. Interpolate the data to have N evenly spaced points along the scaled pseudotime vector using a sliding window of Gaussian weights

2. Determine the genes of interest for alignment

3. Align your trajectory among the selected genes either along the whole trajectory or along a partial segment.

At this point you may either apply the alignment to compare your trajectories or you can use the relative shift in pseudotime to compare differences between your trajectories.

Included in our package is an example: single-cell RNA-seq data set downloaded from GEO(GSE48968) and as shown in figure 2 of the CellAlign publication. They are mouse dendritic cells with two stimulations - LPS and PAM - preprocessed to eliminate bad gene reads and cells with fewer than 2500 genes measured. Here we include a 99 gene module that are induced two-fold by the LPS stimulation for the global alignment and an 89 gene module for the local alignment. We used a technique called diffusion pseudotime [1] to find the trajectories for these data. This pre-processing is typical for single-cell RNA-seq trajectory building, but is by no means the only method of generating such a trajectory. The user should pre-process his or her own data in a technically and biologically relevant way to generate a trajectory for each sample for comparison. The user should keep in mind that the pseudotime metric will be different for each trajectory-finding algorithm because the dimensionality reduction and distance metrics are different for each algorithm [2].  

```{r include = FALSE}
data(expGlobalLPS)
data(expGlobalPAM)
data(expLocalLPS)
data(expLocalPAM)
data(trajLPS)
data(trajPAM)
trajLPS = as.numeric(t(trajLPS))
names(trajLPS) = colnames(expGlobalLPS)
trajPAM = as.numeric(t(trajPAM))
names(trajPAM) = colnames(expGlobalPAM)


```

Interpolation and scaling
-----
The first step is to interpolate the data along the trajectory to represent the data by N (default 200) equally spaced points along the pseudotime trajectory. We included this step because single-cell measurements are often sparse or heterogeneous along the trajectory, leaving gaps that cannot be aligned. Cell-Align interpolates the gene-expression values of equally spaced artificial points using the real single-cell expression data. The expression values of the interpolated points are calculated using all cells, with each single cell assigned a weight given by a Gaussian distribution centered at the interpolated point and a width assigned by a parameter called **winSz**. The default **winSz** is 0.1, as this is the range that preserves the dynamics of the trajectory without including excessive noise for standard single cell data sets. 

```{r fig.height=5}
numPts = 200
interGlobalLPS = cellAlign::interWeights(expDataBatch = expGlobalLPS, trajCond = trajLPS,
                                         winSz = 0.1, numPts = numPts)
interGlobalPAM = cellAlign::interWeights(expDataBatch = expGlobalPAM, trajCond = trajPAM,
                                         winSz = 0.1, numPts = numPts)

```

```{r echo=FALSE}
require(ggplot2)
require(reshape2)
require(pheatmap)
sharedMarkers = intersect(rownames(expGlobalLPS), rownames(expGlobalPAM))
whichgene=sharedMarkers[1]
selectedLPS<-interGlobalLPS$interpolatedVals[whichgene,]
selectedPAM<-interGlobalPAM$interpolatedVals[whichgene,]

dfLPSi = data.frame(traj = interGlobalLPS$traj, value=(selectedLPS), error=interGlobalLPS$error[whichgene,])
dfLPS = data.frame(traj = trajLPS, t(expGlobalLPS[whichgene,]))
dfPAMi = data.frame(traj = interGlobalPAM$traj, value=(selectedPAM), error=interGlobalPAM$error[whichgene,])
dfPAM = data.frame(traj = trajPAM, t(expGlobalPAM[whichgene,]))
dfLPSM = melt(dfLPS, id.vars = 'traj')
dfPAMM = melt(dfPAM, id.vars = 'traj')
#plot of an example gene and its interpolation with error bars
ggplot(dfLPSi, aes(x=traj,y=value)) +  geom_errorbar(aes(ymin=value-error/2, ymax=value+error/2)) + geom_line(size = 2) + geom_point(data=dfLPSM, aes(x=traj,y=value)) + ggtitle(whichgene) 

```

Next, you can scale the expression matrix. This is usual for gene expression analysis, but optional. CellAlign alignment will still run without scaling the expression matrix.

```{r}
#scale the interpolated data (Recommended):
interScaledGlobalLPS = cellAlign::scaleInterpolate(interGlobalLPS)
interScaledGlobalPAM = cellAlign::scaleInterpolate(interGlobalPAM)
```

Choose the genes -
====
Now, the user has to choose which genes to include in the alignment. The choice of genes depends on the goals of the experiment. 

For example, one can compare the expression **gene-by-gene between samples** to check for individual gene expression similarity across the transcriptome. Also cellAlign can compare expression **gene-by-gene within a sample** to identify modules or compare expression dynamics of single genes.

For multiple genes, the most straightforward choice is to use **the genes that were included in the trajectory-finding** preprocessing step, since these genes were likely chosen for their relevance to the biological question of interest. 

In the case of experiment replicates or comparison between different stimulations of the same sample, it might be worthwhile to include all genes. 

# Alignment -

Finally, there is the alignment step. CellAlign operates much like sequence alignment algorithms, quantifying overall similarity in expression throughout the trajectory (global alignment), or finding areas of highly conserved expression (local alignment). The method used is based on dynamic time-warping (DTW) [3, 4] in which a dissimilarity matrix is computed from the distance between each pair of points in the gene-expression space. The user can specify the distance measure to be used - "euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski". Here we use the  Euclidean distance. 

```{r echo=FALSE, message=FALSE, fig.height=5, fig.width=5.5}
A=calcDistMat(interScaledGlobalPAM$scaledData[,1:10],interScaledGlobalLPS$scaledData[,1:10], dist.method = 'Euclidean')
pheatmap(A, cluster_cols = F, cluster_rows=F, main = "LPS vs PAM distances, 1st 10 points",
          show_rownames = F, show_colnames = F,display_numbers = TRUE)

```

Cell-Align then finds a path through the matrix that minimizes the overall distance while adhering to the following constraints: 
 *	for global alignment the alignment must cover the entire extent of both trajectories, always starting in the upper left of the dissimilarity matrix and ending in the lower right.
 
 *	for local alignment the alignment is restricted only to highly similar cells, yielding as output regions with conserved expression dynamics 

Intuitively, the optimal alignment runs along a "valley" within the dissimilarity matrix. 

The cellAlign package performs DTW to align pseudo-temporal sequences by "warping" the temporal axis. DTW is essentially a point-to-point matching algorithm, simliar to DNA sequence alignment algorithms, constrained to maintain the uni-directionality of time and in general without allowing skipped points. At the heart of DTW is a step pattern that assigns penalties to skipping points along the trajectory or excessive repeating of individual points, effectively putting limits on the local slope of the warping curve. CellAlign allows the user to set the step pattern as input (step.pattern = ("symmetric1", "symmetric2", "asymmetric", "rabinerJuang", "symmetric5")). The default pattern, _"symmetric2"_" is adequate for most cellAlign applications because the interpolated data are sufficiently smooth that there is nothing gained from allowing skipped points, and the additional cost for staying on the diagonal balances the fact that the diagonal is the shortest path [5].

**Global Alignment**

The input of the globalAlign function can be either a dissimilarity matrix (as calculated above), or a query and reference matrix as shown below.The pseudotime trajectory vectors must also be included.


```{r warning=FALSE, message=FALSE, fig.height=4, fig.width=4.5}
#perform global alignment of all genes:
alignment = globalAlign(interScaledGlobalPAM$scaledData, interScaledGlobalLPS$scaledData,
                                   scores = list(query = interScaledGlobalPAM$traj, 
                                                 ref = interScaledGlobalLPS$traj),
                                   sigCalc = F, numPerm = 20)
plotAlign(alignment)


```
Here we see the optimal alignment through the dissimilarity matrix. 


** Mapping the results** 

Next we can perform the **mapping** from the interpolation to the real data. This function maps the warping of the interpolated points back onto the original data.
```{r warning=FALSE, message=FALSE, fig.height=4, fig.width=4.5}
mapping = mapRealDataGlobal(alignment,intTrajQuery = interScaledGlobalPAM$traj, realTrajQuery = trajPAM,
                            intTrajRef = interScaledGlobalLPS$traj, realTrajRef = trajLPS)
plotMapping(mapping)

```

**Local Alignment**
The approach for local alignment is different from the global alignment because rather than assessing the overall differences between the reference and query trajectories, cellAlign identifies regions of the trajectory that are closer to each other than a user-provided **similarity threshold**. The threshold gates the dissimilarity matrix, limiting the alignment to local regions with distances smaller than the threshold. The algorithm then backtracks from local minima, stopping when it encounters a point greater than the threshold. It returns the longest aligned trajectory.

```{r warning=FALSE, message=FALSE, fig.height=4, fig.width=4.5}
Thresh=0.2
numPts = 200
interLocalLPS = interWeights(expDataBatch = expLocalLPS, trajCond = trajLPS, winSz = 0.1, numPts = numPts)
interLocalPAM = interWeights(expDataBatch = expLocalPAM, trajCond = trajPAM, winSz = 0.1, numPts = numPts)
interScaledLocalLPS = cellAlign::scaleInterpolate(interLocalLPS)
interScaledLocalPAM = cellAlign::scaleInterpolate(interLocalPAM)

A=calcDistMat(interScaledLocalPAM$scaledData,interScaledLocalLPS$scaledData, dist.method = 'Euclidean')
A[A > 10*Thresh] <- max(A)
alignment = localAlign(interScaledLocalPAM$scaledData,interScaledLocalLPS$scaledData,threshPercent = Thresh)

```
```{r echo=FALSE, message=FALSE, warning=FALSE, fig.height=4,fig.width=4.5}
costMat = t(apply(A,1,function(x){return(as.numeric(x))}))
linearInd = cellAlign::sub2ind(nrow(A), alignment$align[[1]]$index1, alignment$align[[1]]$index2)
costMat[linearInd] = NA
costMat = data.frame(costMat, row.names=1:numPts)
colnames(costMat) = 1:numPts
pheatmap(costMat, cluster_cols = F, cluster_rows=F, border_color = NA,
         main = 'gated search region',
         show_rownames = F, show_colnames = F)

plotAlign(alignment)

BLPS=colMeans(interScaledLocalLPS$scaledData)
BPAM=colMeans(interScaledLocalPAM$scaledData)
p=unique(alignment$align[[1]]$index1)
q=unique(alignment$align[[1]]$index2)
plot(1:200,BPAM,xlab = "pseudotime",ylab = "mean interpolated expression", main = "unaligned mean expression",ylim = c(0,1.1))
points(p,BPAM[p],col="red")
points(1:200,BLPS,col="grey60")
points(q,BLPS[q],col="red")
text(90,1,"LPS")
text(150,1,"PAM")
text(125,.3,"red points are conserved")

```
# Alignment of branched trajectories
Sometimes your data will create a branched trajectory when using trajectory-building algorithms that support this kind of output (such as Monocle2 or diffusion pseudotime). 

# Clustering of aligned trajectories
```{r}
pseudotimeClust <- function(x, y, k = 2)
x - interpolated scaled expression condition 1
y - interpolated scaled expression condition 2
k - number of clusters
```
This function gets interpolated scaled expression of genes (rows) along pseudotime (columns), applies global alignment per gene and cluster the genes by their pseudotime shifts (k means clustering).

# References
1.	Haghverdi, L., et al., Diffusion pseudotime robustly reconstructs lineage branching. Nat Methods, 2016. 13(10): p. 845-8.
2.	Cannoodt, R., W. Saelens, and Y. Saeys, Computational methods for trajectory inference from single-cell transcriptomics. Eur J Immunol, 2016. 46(11): p. 2496-2506.
3.	Aach, J. and G.M. Church, Aligning gene expression time series with time warping algorithms. Bioinformatics, 2001. 17(6): p. 495-508.
4.	Sankoff, D. and J.B. Kruskal, Time warps, string edits, and macromolecules : the theory and practice of sequence comparison. 1983, Reading, Mass.: Addison-Wesley Pub. Co., Advanced Book Program. xii, 382 p.
5. Jiaping Zhao, Laurent Itti, shapeDTW: shape Dynamic Time Warping. 2016, arXiv:1606.01601 [cs.CV]


