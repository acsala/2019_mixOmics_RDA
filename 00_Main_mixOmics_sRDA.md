
# Table of Contents

1.  [Summary](#org63071ca)
2.  [Introduction to sRDA in mixOmics](#org8fc0fe9)
3.  [Comparison in terms of objective functions](#org94cdf7e)
4.  [Compare sRDA and sPLS canonical mode with 5 components on mixOmics breast.TCGA data](#org360ae3c)
    1.  [Overlapping variables](#orgc51c8a7)



<a id="org63071ca"></a>

# Summary

Redundancy analysis (RDA) is a latent variable based regression
method, often referred as the multivariate equivalent of multiple
regression (Legendre 2012). RDA was first described by Wollenberg
(van den Wollenberg 1977), its development was motivated by
expending the ordinary least squares multiple regression to multiple
response variables. It has been extensively used in ecology and
chemo-metrics, and lately has been developed for high-dimensional
data analysis in the form of sparse Redundancy Analysis (sRDA).

sRDA can be used when the objective is to find a combination of
predictor variables (a latent variable) that explain the most
variance in the response variables. In omics data analyses, for
example, one might be interested to find a combination of multiple
mRNA levels that explain the most variation in multiple proteomics
levels, measured in the sample of patients with a particular
disease. After applying sRDA in this setting, the extracted
variables can be interpreted as bio-markers for the given
disease. Moreover, the explanatory-response properties of the
variables are preserved, thus indicating that changing (some of)
those mRNA levels that form the latent variable would result in a
change in the proteomics levels. sRDA also provides an indication of
the strength of "dependency" of all response variables on the latent
variable in terms of correlation coefficients, which can be
interpreted as the multiple regression correlation coefficient per
response variable.

In this document, we describe when and how to use sRDA for omics
data analysis. We describe the differences between sRDA and other
multivariate methods available from mixOmics in terms of objective
function and also in terms of bilogical interpretability.


<a id="org8fc0fe9"></a>

# Introduction to sRDA in mixOmics

sRDA can be applied to answer the following biological
questions:

I have two omics data sets measured on the same patients, one data
set measured on mRNA levels and a second data set measured on
protein levels. I assume that mRNAs are predictors for protein
levels and I am interested in which are the "important" mRNAs that
explain the most change in the protein levels? Which are the
proteins that response most and therefore are most affected by these
mRNAs? What is the quantified effect of the "important" mRNAs?

Following this example, we demonstrate the application of sRDA to
a filtered breast cancer data from he Cancer Genome Atlas (TCGA,
<http://cancergenome.nih.gov/>) available in mixOmics. The expression
levels of 200 mRNA and the abundance of 142 proteins are available
on 150 patients with breast cancer. 

    # load cancer data from mixOmics package
    data(breast.TCGA)
    
    # extract scaled training data (with mean = 0 and 
    # standard deviation = 1 columns)
    X <- scale(breast.TCGA$data.train$mrna)
    Y <- scale(breast.TCGA$data.train$protein)
    
    #call the sRDA function to run sRDA on data sets X and Y
    res_sRDA <- get_PLS_CCA_RDA_results(X = X, Y = Y,
    				    nr_nonz = 20, 
    				    nr_comp = 2,
    				    penalty_mode = "enet",
    				    CCA = F)$res_sRDA
    data.frame(dim(X),dim(Y))

<table border="2" cellspacing="0" cellpadding="6" rules="groups" frame="hsides">


<colgroup>
<col  class="org-right" />

<col  class="org-right" />
</colgroup>
<thead>
<tr>
<th scope="col" class="org-right">dim.X.</th>
<th scope="col" class="org-right">dim.Y.</th>
</tr>
</thead>

<tbody>
<tr>
<td class="org-right">150</td>
<td class="org-right">150</td>
</tr>


<tr>
<td class="org-right">200</td>
<td class="org-right">142</td>
</tr>
</tbody>
</table>

<table border="2" cellspacing="0" cellpadding="6" rules="groups" frame="hsides">


<colgroup>
<col  class="org-left" />

<col  class="org-right" />

<col  class="org-left" />

<col  class="org-right" />
</colgroup>
<thead>
<tr>
<th scope="col" class="org-left">mRNA<sub>results</sub><sub>names</sub></th>
<th scope="col" class="org-right">&#xa0;</th>
<th scope="col" class="org-left">protein<sub>results</sub><sub>names</sub></th>
<th scope="col" class="org-right">&#xa0;</th>
</tr>
</thead>

<tbody>
<tr>
<td class="org-left">CCNA2</td>
<td class="org-right">0.152</td>
<td class="org-left">ER-alpha</td>
<td class="org-right">-0.813</td>
</tr>


<tr>
<td class="org-left">ZNF552</td>
<td class="org-right">-0.126</td>
<td class="org-left">GATA3</td>
<td class="org-right">-0.768</td>
</tr>


<tr>
<td class="org-left">KDM4B</td>
<td class="org-right">-0.107</td>
<td class="org-left">Cyclin<sub>B1</sub></td>
<td class="org-right">0.745</td>
</tr>


<tr>
<td class="org-left">PREX1</td>
<td class="org-right">-0.107</td>
<td class="org-left">ASNS</td>
<td class="org-right">0.716</td>
</tr>


<tr>
<td class="org-left">LRIG1</td>
<td class="org-right">-0.093</td>
<td class="org-left">PR</td>
<td class="org-right">-0.704</td>
</tr>


<tr>
<td class="org-left">E2F1</td>
<td class="org-right">0.052</td>
<td class="org-left">Cyclin<sub>E1</sub></td>
<td class="org-right">0.704</td>
</tr>


<tr>
<td class="org-left">C4orf34</td>
<td class="org-right">-0.05</td>
<td class="org-left">AR</td>
<td class="org-right">-0.692</td>
</tr>


<tr>
<td class="org-left">ASPM</td>
<td class="org-right">0.047</td>
<td class="org-left">JNK2</td>
<td class="org-right">-0.654</td>
</tr>


<tr>
<td class="org-left">NTN4</td>
<td class="org-right">-0.038</td>
<td class="org-left">INPP4B</td>
<td class="org-right">-0.654</td>
</tr>


<tr>
<td class="org-left">FUT8</td>
<td class="org-right">-0.033</td>
<td class="org-left">CDK1</td>
<td class="org-right">0.65</td>
</tr>


<tr>
<td class="org-left">TTC39A</td>
<td class="org-right">-0.03</td>
<td class="org-left">Bcl-2</td>
<td class="org-right">-0.558</td>
</tr>


<tr>
<td class="org-left">STC2</td>
<td class="org-right">-0.025</td>
<td class="org-left">ER-alpha<sub>pS118</sub></td>
<td class="org-right">-0.527</td>
</tr>


<tr>
<td class="org-left">SLC19A2</td>
<td class="org-right">-0.022</td>
<td class="org-left">PDK1<sub>pS241</sub></td>
<td class="org-right">-0.515</td>
</tr>


<tr>
<td class="org-left">MEX3A</td>
<td class="org-right">0.017</td>
<td class="org-left">p53</td>
<td class="org-right">0.512</td>
</tr>


<tr>
<td class="org-left">C18orf1</td>
<td class="org-right">-0.016</td>
<td class="org-left">Chk2</td>
<td class="org-right">0.51</td>
</tr>


<tr>
<td class="org-left">MTL5</td>
<td class="org-right">-0.015</td>
<td class="org-left">DJ-1</td>
<td class="org-right">-0.5</td>
</tr>


<tr>
<td class="org-left">NCAPG2</td>
<td class="org-right">0.014</td>
<td class="org-left">Chk2<sub>pT68</sub></td>
<td class="org-right">0.475</td>
</tr>


<tr>
<td class="org-left">MED13L</td>
<td class="org-right">-0.01</td>
<td class="org-left">Caveolin-1</td>
<td class="org-right">-0.466</td>
</tr>


<tr>
<td class="org-left">FAM63A</td>
<td class="org-right">-0.009</td>
<td class="org-left">P-Cadherin</td>
<td class="org-right">0.463</td>
</tr>


<tr>
<td class="org-left">SEMA3C</td>
<td class="org-right">-0.006</td>
<td class="org-left">p27<sub>pT198</sub></td>
<td class="org-right">0.46</td>
</tr>
</tbody>
</table>

In this example, we applied sRDA to extract two latent variables
from the mRNA measurements (data set X) to explain the variance in
the protein levels (data set Y). Note that we scale our data sets
(i.e. column transformation of mean = 0 and standard deviation = 1).

We explicitly defined that we wish to obtain 20 mRNAs that best
explain the variance in all the protein variables. In order to
analyse the results, we can print out the loadings for the mRNAs and
proteins, which indicates the strength of the particular mRNAs in
the latent variable, and gives the multiple regression coefficient
(in respect to the latent variable) for each protein variables (we
printed out the 20 highest here).

    ER_lm <- lm(Y[,"ER-alpha"] ~ res_sRDA$variates$X[,1])
        plot(res_sRDA$variates$X[,1], Y[,"ER-alpha"], 
    	 ylab = "ER-alpha level", xlab = "Latent variable", 
    	 pch = 16, cex = 1.3)
    abline(ER_lm,col = "red", lwd = 2.5)

<./plots/plot_regression_interpretation_sRDA.pdf>

We can read from the results that ER-alpha protein levels have a
strong negative correlation with the latent variable (with
correlation coefficient -0.813), thus one unit change in the latent
variable correlates with a .813 decrees in the (scaled) protein
levels. We can also read that mRNA CCNA2 has a 0.152 coefficient
with the latent variable, thus a unit change in CCNA2 correlates
with a .152 increase in the latent variable.

Additionally, we can use the mixOmics' plots, such as
plotLoadings(), plotIndiv(), plotVar() or cim() (the detailed
descriptions for these plots can be found in mixOmics' vignette).

    
    plotLoadings(res_sRDA, comp = 1, size.name = rel(0.5))

<./plots/plot_Loadings_sRDA.pdf>

    plotIndiv(res_sRDA)

<./plots/plot_Indiv_sRDA.pdf>

    
    plotVar(res_sRDA)

<./plots/plot_var_sRDA.pdf>

    
    cim(res_sRDA, comp = 1)

<./plots/plot_cim_sRDA.pdf>


<a id="org94cdf7e"></a>

# Comparison in terms of objective functions

Partial least squares (PLS), Canonical correlation analysis (CCA),
Principal component analysis (PCA) and RDA are closely related
multivariate latent variable based methods. In this section, we
describe their relationship in terms of objective functions, which
also applies to the penalized and sparse versions of these methods.

Given \(\mathbf X \in \mathbb{R}^{n \times p}\), a centered and scaled
(with columns of mean zero and standard deviation one, respectively)
predictor data set and \(\mathbf Y \in \mathbb{R}^{n \times q}\), a
centered and scaled response data set, and if our goal is to look
for the first pair of loading vectors, \(\mathbf w \in \mathbb R^p\)
for \(\mathbf X\) and \(\mathbf v \in \mathbb R^q\) for \(\mathbf Y\),
then these methods maximize the following quantities:

\begin{align}
\mathrm{PCA:}&\quad \operatorname{Var}(\mathbf{Xw}) \\
\mathrm{CCA:}&\quad \phantom{\operatorname{Var}(\mathbf {Xw})\cdot {}}
\operatorname{Corr}^2(\mathbf {Xw},\mathbf {Yv})\\
\mathrm{RDA:}&\quad \phantom{\operatorname{Var}(\mathbf {Xw})\cdot{}}
\operatorname{Corr}^2(\mathbf{Xw},\mathbf {Yv})\cdot\operatorname{Var}(\mathbf{Yv})  \\
\mathrm{PLS:}&\quad \operatorname{Var}(\mathbf{Xw})\cdot\operatorname{Corr}^2(\mathbf{Xw},
\mathbf {Yv})\cdot\operatorname{Var}(\mathbf {Yv}) = \operatorname{Cov}^2(\mathbf{Xw},\mathbf {Yv})
\end{align}

PCA is a unsupervised method which aims to maximize the variance of
the linear combinations of all variables from one data set (Eq. (1)). These
linear combinations are also called latent variables. By maximising
the variance of latent variables, the aim is to capture all the
variation in the original data set in the latent variables. By doing
so, one aims to preserve or extract all information from the
original data in a lower dimensional space.

CCA and PLS are close related to PCA, they aim to extract latent
variables that explain variation in the original data. They are
applied to two (or multiple) data sets. CCA aims to maximize the
squared correlations between the latent variables (Eq. (2)), which
results of a linear combination of variables from **X** that has the
maximum correlation with a linear combination of variables from
**Y**. PLS aims to maximize the squared covariance between linear the
latent variables (Eq. (4)), which is a scaled version of CCA's
objective function. Both PLS and CCA are non-directional methods,
the optimum of the objective function doesn't change by
interchanging the input data sets. These methods can be used to
explore linear combinations of variables that have the highest
correlation or covariance (respectively) with each other. CCA is
preferred if the relationship wished to be characterized in
correlations instead of covariances (one is standardized the other
is in the scale of the original variables).

RDA is similar to the aforementioned methods, in the sense that it is
based on latent variable extraction. The objective function of RDA
is to maximise the squared correlation between a linear combination
of the predictor variables and the response data set. Therefore, RDA
is a directional method that distinguishes between the predictor and
response data sets and the optimum of its objective function changes
if one interchanges the input data sets. RDA and can be used to
explore linear combinations of explanatory variables that explain
the most variance in an outcome data set.


<a id="org360ae3c"></a>

# Compare sRDA and sPLS canonical mode with 5 components on mixOmics breast.TCGA data

In order to compare sRDA to sPLS, Both methods are applied on a
dataset with 150 patients and 200 mrna and 142 protein measurement.

5 components are extracted with enforcing 10 non-zero variables (no
optimization procedure).


<a id="orgc51c8a7"></a>

## Overlapping variables

    set.seed(88)
    ncomp <- 5
    nr_nonz <- 10
    res_all <- get_PLS_CCA_RDA_results(X = X, Y = Y, 
    				 nr_nonz = nr_nonz, 
    				 nr_comp = ncomp,
    				 pls_mode = "canonical",
    				 penalty_mode = "ust",
    				 CCA = F)
    
    res_sRDA <- res_all$res_sRDA
    res_spls <- res_all$res_spls
    
    # Compare the variables selection
    res_sRDA <- get_nonzero_variables(res_object = res_sRDA)
    res_spls <- get_nonzero_variables(res_object = res_spls)
    
    #res_sRDA$nz_loading_names[["X"]]
    #res_spls$nz_loading_names[["X"]]
    
    plot(1:ncomp,get_nr_of_common_components(res_sRDA, res_spls, 
    				       ncomp = ncomp), 
       ylab = "Nr of common variables", xlab = "Component", 
       ylim = c(0,nr_nonz), pch = 19, lwd = 2)

<./plots/fig_overlap.pdf>

The variables selected in the first component
are identical for PLS and RDA, there are 5 overlapping variables in
the second component, and there are no overlapping variables in
consecutive latent variables.

