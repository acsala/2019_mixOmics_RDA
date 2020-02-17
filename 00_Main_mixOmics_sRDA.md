
# Table of Contents

1.  [Introduction to sRDA in mixOmics](#org7a4e0f7)
    1.  [Biological questions to answer with sRDA](#orge43e88c)
        1.  [Which are the "important" explanatory variables that explain the most variance in the response variables?](#org130606e)
        2.  [Which response variables are affected most strongly by the "important" explanatory variables?](#org6929597)
        3.  [What is the "quantified" effect of the explanatory variables on the response variables](#org119c9ef)
2.  [Comparison of latent variable based methods in terms of objective functions](#org1a2b0bc)
3.  [Cross validated analysis](#org1279ee3)
4.  [Compare sRDA and sPLS canonical mode with 5 components on mixOmics breast.TCGA data](#org5839204)
    1.  [Overlapping variables](#orgfd491ff)
    2.  [Benchmark correlation](#orgea94e63)
5.  [References](#org6fd3c3b)



<a id="org7a4e0f7"></a>

# Introduction to sRDA in mixOmics

Redundancy analysis (RDA) is a directional, latent variable based
regression method, often referred to as the multivariate equivalent
of multiple regression <sup id="59e1f790b3b37a579e8a75fcebe25b6a"><a href="#legendre12_numer" title="Legendre, Numerical ecology, Elsevier (2012).">legendre12_numer</a></sup>. RDA accounts for the
dependency between data sets, which distinguishes it from other
latent variable based methods in the mixOmics package
<sup id="96148de61680d5941b2957b7ed478649"><a href="#Rohart_2017" title="Rohart, Gautier, Singh, \&amp; L&#234; Cao, mixOmics: an R package for &#8216;omics feature selection  and multiple data integration, v(), (2017).">Rohart_2017</a></sup>. RDA can be used when the goal of the
analysis is to find a combination of \(\color{red}explanatory\
  variables\) from a explanatory data set that explain the most
variance in all the \(\color{blue}response\ variables\) in a response
data set (see Figure [fig:conceputal_framework_RDA](#fig:conceputal_framework_RDA)).

<./plots/plot_sRDA_concept.pdf>

RDA was first described by Wollenberg
<sup id="789070d7c6e50336d15a21f274c89c09"><a href="#wollenberg77_redun_analy_alter_canon_correl_analy" title="Arnold van den Wollenberg, Redundancy Analysis an Alternative for Canonical  Correlation Analysis, {Psychometrika}, v(2), 207-219 (1977).">wollenberg77_redun_analy_alter_canon_correl_analy</a></sup>, its
development was motivated by expending the ordinary least squares
multiple regression to multiple response variables. It has been
extensively used in ecology and chemo-metrics, and lately has been
developed for high-dimensional data analysis in the form of sparse
Redundancy Analysis (sRDA)
<sup id="294ead7176458d8a872380daa50f5068"><a href="#csala17_spars_redun_analy_high_dimen" title="Attila Csala, Frans P J M Voorbraak, Aeilko H, Zwinderman \&amp; Michel H Hof, Sparse Redundancy Analysis of High-Dimensional  Genetic and Genomic Data, {Bioinformatics}, v(20), 3228-3234 (2017).">csala17_spars_redun_analy_high_dimen</a></sup>. sRDA can be used when the
objective is to find a combination of the subset of
\(\color{red}explanatory\ variables\) (also called the latent
variable) that explain the most variance in the
\(\color{blue}response\ variables\). In omics data analyses, for
example, one might be interested to find a combination of
\(\color{red}multiple\ mRNA\ levels\) that explain the most variation
in \(\color{blue}multiple\ proteomics\ levels\), measured in the
sample of patients with a particular disease. Applying sRDA in this
setting, the extracted variables can be interpreted as bio-markers
for the given disease. Moreover, the
\(\color{red}explanatory-\color{blue}response\) properties of the
variables are preserved, thus indicating that changing (some of)
the mRNA levels that form the latent variable would result in a
change in the proteomics levels. sRDA also provides an indication of
the strength of "dependency" of all response variables on the latent
variable in terms of correlation coefficients, which can be
interpreted as the multiple regression correlation coefficient per
response variable.


<a id="orge43e88c"></a>

## Biological questions to answer with sRDA

sRDA can be applied to answer the following biological
questions:

I have two omics data sets measured on the same patients, one data
set measured on mRNA levels and a second data set measured on
protein levels. 

1.  I assume that mRNAs are \(\color{red}explanatory\ variables\) for
    protein levels (i.e. \(\color{blue}response\ variables\)). I am
    interested to find out which are the "important" mRNAs that
    explain the most variance in the protein levels?
2.  And which are the proteins that response most strongly and
    therefore are most affected by the "important" mRNAs?
3.  And what is the quantified effect of the "important" mRNAs on the
    proteins?

As a motivating example, we apply sRDA to a filtered breast cancer
data from he Cancer Genome Atlas (TCGA,
<http://cancergenome.nih.gov/>) to answer these questions. The
filtered TCGA data is available in mixOmics, here are
\(\color{red}expression\ levels\ of\ 200\ mRNA\) and the
\(\color{blue}abundance\ of\ 142\ proteins\) are available on 150
patients with breast cancer.

    #load sRDA and helper functions for mixOmics from Github
    library(sRDA)
    library(devtools)
    url <- "https://raw.githubusercontent.com/acsala/2019_mixOmics_RDA/master/00_helper_functions.R"
    source_url(url)

    # load cancer data from mixOmics package
    library(mixOmics)
    data(breast.TCGA)
    
    #set seed for reproducibility
    set.seed(100)
    
    X <- breast.TCGA$data.train$mrna
    Y <- breast.TCGA$data.train$protein
    
    #call the sRDA function to run sRDA on data sets X and Y
    res_sRDA <- sRDA_mixOmics(X = X, Y = Y,
    			  keepX = 20,
    			  penalty_mode = "enet",
    			  ncomp = 5,
    			  scale = TRUE)

First we load sRDA from CRAN and its mixOMICS helper functions from
Github (Listing [code_load_RDA_helpers](#code_load_RDA_helpers)) and then apply sRDA to
the TCGA data from mixOmics (Listing [code_run_sRDA_on_TCGA](#code_run_sRDA_on_TCGA)).

In this example, we applied sRDA to extract two latent variables
(LVs) from the mRNA measurements (data set X) to explain the
variance in the protein levels (data set Y). The two LVs represent
independent sets of the combination of 20 mRNAs that explain the
most variance in the protein levels. With these results, we can
answer our biological questions.


<a id="org130606e"></a>

### Which are the "important" explanatory variables that explain the most variance in the response variables?

We explicitly defined that we wish to obtain 20 mRNAs that best
explain the variance in all the protein variables (Listing
[code_run_sRDA_on_TCGA](#code_run_sRDA_on_TCGA)). In order to analyse the results, we
can print out the loadings for the mRNAs and proteins, which
indicates the strength of the particular mRNAs in the latent
variable, and gives the multiple regression coefficient (in
respect to the latent variable) for each protein variables (see
Table [tab:sRDA_first_LV](#tab:sRDA_first_LV) for the selected variables and their
loadings, where we printed out only the top 20 absolute loadings
of the protein variables).

1.  I assume that mRNAs are \(\color{red}explanatory\ variables\) for
    protein levels (i.e. \(\color{blue}response\ variables\)). I am
    interested to find out which are the "important" mRNAs that
    explain the most change in the protein levels?
    
    We can see these important mRNAs in Table [tab:sRDA_first_LV](#tab:sRDA_first_LV).

    
    lv = 1
    protein_names <-
      names(sort(abs(res_sRDA$loadings$Y[,lv]), decreasing = T)[1:20])
    mRNA_names <-
      names(sort(abs(res_sRDA$loadings$X[,lv]), decreasing = T)[1:20])
    
    print(cbind(mRNA_names,
    	    round(res_sRDA$loadings$X[,lv][mRNA_names],
    		  digits = 3), 
    	    protein_names,
    	    round(res_sRDA$loadings$Y[,lv][protein_names],
    		  digits = 3)))

<table border="2" cellspacing="0" cellpadding="6" rules="groups" frame="hsides">
<caption class="t-above"><span class="table-number">Table 1:</span> mRNA and protein variables extracted in the first latent variable from the TCGA data set with sRDA <a name="tab:sRDA_first_LV"></a></caption>

<colgroup>
<col  class="org-left" />

<col  class="org-right" />

<col  class="org-left" />

<col  class="org-right" />
</colgroup>
<thead>
<tr>
<th scope="col" class="org-left">mRNA<sub>names</sub></th>
<th scope="col" class="org-right">&#xa0;</th>
<th scope="col" class="org-left">protein<sub>names</sub></th>
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


<a id="org6929597"></a>

### Which response variables are affected most strongly by the "important" explanatory variables?

1.  And which are the proteins that response most strongly and
    therefore are most affected by the "important" mRNAs?
    
    We can find these proteins by ordering our results based on the
    absolute loadings of the proteins. The top 20 most affected
    proteins can be found in Table [tab:sRDA_first_LV](#tab:sRDA_first_LV).


<a id="org119c9ef"></a>

### What is the "quantified" effect of the explanatory variables on the response variables

We can find regression weights for the mRNA variables (i.e. the
regression weight between the mRNA and the latent variable) and
the) multiple regression correlation coefficient for the protein
levels (i.e. the multiple regression coefficient between the
protein levels and the LV) in Table [tab:sRDA_first_LV](#tab:sRDA_first_LV) that
helps us to quantify the relationships we found.  ER-alpha protein
levels have a strong negative correlation with the latent variable
(with correlation coefficient -0.813), thus one unit change in the
latent variable correlates with a .813 decrees in the (scaled)
protein levels (see Figure
[fig_regression_interpretation_sRDA](#fig_regression_interpretation_sRDA)). We can also read that
mRNA CCNA2 has a 0.152 coefficient with the latent variable, thus
a unit change in CCNA2 correlates with a .152 increase in the
latent variable (see Table [tab:sRDA_first_LV](#tab:sRDA_first_LV)).

    ER_lm <- lm(Y[,"ER-alpha"] ~ res_sRDA$variates$X[,1])
    
    plot_cus(res_sRDA$variates$X[,1], Y[,"ER-alpha"], 
    	 ylab = "ER-alpha level", xlab = "Latent variable",
    	 ylim = c(-4,4),
    	 bg = col[1])
    
    abline(c(ER_lm), col=c(regression_line_col), lwd=c(2.5))

<./plots/fig_regression_interpretation_sRDA.pdf>

1.  And what is the quantified effect of the "important" mRNAs on the
    proteins?
    
    We can read from the results that ER-alpha protein levels have a
    strong negative correlation with the latent variable (with
    correlation coefficient -0.813), thus one unit change in the
    latent variable correlates with a .813 decrees in the (scaled)
    protein levels (see Figure
    [fig_regression_interpretation_sRDA](#fig_regression_interpretation_sRDA)). We can also read that
    mRNA CCNA2 has a 0.152 coefficient with the latent variable,
    thus a unit change in CCNA2 correlates with a .152 increase in
    the latent variable (see Table [tab:sRDA_first_LV](#tab:sRDA_first_LV)).

Additionally, we can use the mixOmics' plots, such as
plotLoadings(), plotIndiv(), plotVar() or cim() (the detailed
descriptions for these plots can be found in mixOmics' vignette).

    
    plotLoadings(res_sRDA, ndisplay = 20, comp = 1, size.name = rel(0.9),
    	     contrib = "max", size.title = rel(1.3))

<./plots/plot_Loadings_sRDA.pdf>

    plotIndiv(res_sRDA)

<./plots/plot_Indiv_sRDA.pdf>

    
    plotVar(res_sRDA, cutoff = 0.5)

<./plots/plot_var_sRDA.pdf>

    cim(res_sRDA, scale = T, threshold = 0.35)

<./plots/plot_cim_sRDA.pdf>

    network(res_sRDA,
    	color.node = c("orange","blue"),
    	cutoff = 0.6,
    	row.names = T, 
    	col.names = T,
    	comp = 1)

<./plots/network_sRDA.pdf>

In this section we described when and how to use sRDA for omics data
analysis. In the following we describe the differences between sRDA
and other multivariate methods available from mixOmics in terms of
objective function and also in terms of bilogical interpretability.


<a id="org1a2b0bc"></a>

# Comparison of latent variable based methods in terms of objective functions

Partial least squares (PLS), Canonical correlation analysis (CCA),
Principal component analysis (PCA) and RDA are closely related
multivariate latent variable based methods, all available in the
mixOmics package. In this section, we
describe their relationship in terms of objective functions, which
also applies to the penalized and sparse versions of these methods.

Given \(\mathbf X \in \mathbb{R}^{n \times p}\), a centered and scaled
(with columns of mean zero and standard deviation one, respectively)
predictor data set and \(\mathbf Y \in \mathbb{R}^{n \times q}\), a
centered and scaled response data set. Our goal is to look for the
first pair of loading vectors, \(\mathbf w \in \mathbb R^p\) for
\(\mathbf X\) and \(\mathbf v \in \mathbb R^q\) for \(\mathbf Y\), then
these methods maximize the following quantities:

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
the linear combinations of all variables from one data set
(Eq. (1)). These linear combinations are also called latent
variables. By maximising the variance of latent variables, the aim
is to capture all the variation in the original data set in the
latent variables. By doing so, often nearly all information from the
original data can be preserved in a lower dimensional latent
variable space.

CCA and PLS are close related to PCA, they aim to extract latent
variables that explain variation in the original data. They are
applied to two (or multiple) data sets. CCA aims to maximize the
squared correlations between the latent variables (Eq. (2)), which
results of a linear combination of variables from **X** that has the
maximum correlation with a linear combination of variables from
**Y**. PLS aims to maximize the squared covariance between linear the
latent variables (Eq. (4)). Both PLS and CCA are non-directional
methods, that is the optimum of the objective function doesn't
change by interchanging the input data sets. These methods can be
used to explore linear combinations of variables that have the
highest correlation or covariance (respectively) with each
other. CCA is preferred if the relationship wished to be
characterized in correlations instead of covariances (one is
standardized the other is in the scale of the original variables).

RDA is similar to the aforementioned methods, in the sense that it
is based on latent variable extraction. The objective function of
RDA is to maximise the squared correlation between a linear
combination of the predictor variables and the response data set
(Eq. (3)). Therefore, RDA is a directional method that distinguishes
between the predictor and response data sets and the optimum of its
objective function changes if one interchanges the input data
sets. RDA and can be used to explore linear combinations of
explanatory variables that explain the most variance in an outcome
data set.


<a id="org1279ee3"></a>

# Cross validated analysis

We can (and we should) cross-validate for the optimum penalty
parameters. There is Elastic net (enet) and Univariate soft
thresholding (ust) available for sRDA
(<sup id="294ead7176458d8a872380daa50f5068"><a href="#csala17_spars_redun_analy_high_dimen" title="Attila Csala, Frans P J M Voorbraak, Aeilko H, Zwinderman \&amp; Michel H Hof, Sparse Redundancy Analysis of High-Dimensional  Genetic and Genomic Data, {Bioinformatics}, v(20), 3228-3234 (2017).">csala17_spars_redun_analy_high_dimen</a></sup>), below we use enet with
a Ridge penalty of 0.1 and a grid of non-zeros between 5 and 75 (the
longer the grid the longer it takes to run the cross validation)
(Listing [code_sRDA_cross_validate](#code_sRDA_cross_validate)). Under this setting the
cross-validation ran for 2.5 minutes and selected 25 as the optimum
non-zero parameter for the mRNA data source (Figure
[fig_sRDA_cross_validation](#fig_sRDA_cross_validation)).

<table border="2" cellspacing="0" cellpadding="6" rules="groups" frame="hsides">
<caption class="t-above"><span class="table-number">Table 2:</span> sRDA code and results for cross validated analysis <a name="code_sRDA_cross_validate"></a></caption>

<colgroup>
<col  class="org-left" />

<col  class="org-left" />

<col  class="org-right" />

<col  class="org-left" />

<col  class="org-right" />
</colgroup>
<thead>
<tr>
<th scope="col" class="org-left">&#xa0;</th>
<th scope="col" class="org-left">mRNA<sub>names</sub></th>
<th scope="col" class="org-right">&#xa0;</th>
<th scope="col" class="org-left">protein<sub>names</sub></th>
<th scope="col" class="org-right">&#xa0;</th>
</tr>
</thead>

<tbody>
<tr>
<td class="org-left">CCNA2</td>
<td class="org-left">CCNA2</td>
<td class="org-right">0.15</td>
<td class="org-left">ER-alpha</td>
<td class="org-right">-0.817</td>
</tr>


<tr>
<td class="org-left">ZNF552</td>
<td class="org-left">ZNF552</td>
<td class="org-right">-0.123</td>
<td class="org-left">GATA3</td>
<td class="org-right">-0.772</td>
</tr>


<tr>
<td class="org-left">PREX1</td>
<td class="org-left">PREX1</td>
<td class="org-right">-0.108</td>
<td class="org-left">Cyclin<sub>B1</sub></td>
<td class="org-right">0.744</td>
</tr>


<tr>
<td class="org-left">KDM4B</td>
<td class="org-left">KDM4B</td>
<td class="org-right">-0.103</td>
<td class="org-left">ASNS</td>
<td class="org-right">0.725</td>
</tr>


<tr>
<td class="org-left">LRIG1</td>
<td class="org-left">LRIG1</td>
<td class="org-right">-0.095</td>
<td class="org-left">Cyclin<sub>E1</sub></td>
<td class="org-right">0.705</td>
</tr>


<tr>
<td class="org-left">E2F1</td>
<td class="org-left">E2F1</td>
<td class="org-right">0.063</td>
<td class="org-left">PR</td>
<td class="org-right">-0.703</td>
</tr>


<tr>
<td class="org-left">ASPM</td>
<td class="org-left">ASPM</td>
<td class="org-right">0.057</td>
<td class="org-left">AR</td>
<td class="org-right">-0.687</td>
</tr>


<tr>
<td class="org-left">NTN4</td>
<td class="org-left">NTN4</td>
<td class="org-right">-0.057</td>
<td class="org-left">INPP4B</td>
<td class="org-right">-0.664</td>
</tr>


<tr>
<td class="org-left">C4orf34</td>
<td class="org-left">C4orf34</td>
<td class="org-right">-0.053</td>
<td class="org-left">JNK2</td>
<td class="org-right">-0.657</td>
</tr>


<tr>
<td class="org-left">SLC19A2</td>
<td class="org-left">SLC19A2</td>
<td class="org-right">-0.04</td>
<td class="org-left">CDK1</td>
<td class="org-right">0.653</td>
</tr>


<tr>
<td class="org-left">STC2</td>
<td class="org-left">STC2</td>
<td class="org-right">-0.037</td>
<td class="org-left">Bcl-2</td>
<td class="org-right">-0.568</td>
</tr>


<tr>
<td class="org-left">MTL5</td>
<td class="org-left">MTL5</td>
<td class="org-right">-0.037</td>
<td class="org-left">ER-alpha<sub>pS118</sub></td>
<td class="org-right">-0.522</td>
</tr>


<tr>
<td class="org-left">FUT8</td>
<td class="org-left">FUT8</td>
<td class="org-right">-0.037</td>
<td class="org-left">p53</td>
<td class="org-right">0.519</td>
</tr>


<tr>
<td class="org-left">C18orf1</td>
<td class="org-left">C18orf1</td>
<td class="org-right">-0.037</td>
<td class="org-left">Chk2</td>
<td class="org-right">0.515</td>
</tr>


<tr>
<td class="org-left">NCAPG2</td>
<td class="org-left">NCAPG2</td>
<td class="org-right">0.035</td>
<td class="org-left">PDK1<sub>pS241</sub></td>
<td class="org-right">-0.512</td>
</tr>


<tr>
<td class="org-left">MEX3A</td>
<td class="org-left">MEX3A</td>
<td class="org-right">0.032</td>
<td class="org-left">DJ-1</td>
<td class="org-right">-0.495</td>
</tr>


<tr>
<td class="org-left">TTC39A</td>
<td class="org-left">TTC39A</td>
<td class="org-right">-0.032</td>
<td class="org-left">Chk2<sub>pT68</sub></td>
<td class="org-right">0.488</td>
</tr>


<tr>
<td class="org-left">MED13L</td>
<td class="org-left">MED13L</td>
<td class="org-right">-0.028</td>
<td class="org-left">P-Cadherin</td>
<td class="org-right">0.478</td>
</tr>


<tr>
<td class="org-left">FAM63A</td>
<td class="org-left">FAM63A</td>
<td class="org-right">-0.024</td>
<td class="org-left">Caveolin-1</td>
<td class="org-right">-0.476</td>
</tr>


<tr>
<td class="org-left">SNED1</td>
<td class="org-left">SNED1</td>
<td class="org-right">-0.023</td>
<td class="org-left">p27<sub>pT198</sub></td>
<td class="org-right">0.469</td>
</tr>


<tr>
<td class="org-left">SEMA3C</td>
<td class="org-left">SEMA3C</td>
<td class="org-right">-0.017</td>
<td class="org-left">Caspase-7<sub>cleavedD198</sub></td>
<td class="org-right">0.469</td>
</tr>


<tr>
<td class="org-left">PLCD4</td>
<td class="org-left">PLCD4</td>
<td class="org-right">-0.012</td>
<td class="org-left">Smad3</td>
<td class="org-right">-0.464</td>
</tr>


<tr>
<td class="org-left">HTRA1</td>
<td class="org-left">HTRA1</td>
<td class="org-right">-0.002</td>
<td class="org-left">S6</td>
<td class="org-right">0.435</td>
</tr>


<tr>
<td class="org-left">ZEB1</td>
<td class="org-left">ZEB1</td>
<td class="org-right">-0.002</td>
<td class="org-left">Notch1</td>
<td class="org-right">0.432</td>
</tr>


<tr>
<td class="org-left">SLC43A3</td>
<td class="org-left">SLC43A3</td>
<td class="org-right">0</td>
<td class="org-left">Syk</td>
<td class="org-right">0.432</td>
</tr>
</tbody>
</table>

<./plots/fig_sRDA_cross_validation.pdf>


<a id="org5839204"></a>

# Compare sRDA and sPLS canonical mode with 5 components on mixOmics breast.TCGA data

In order to compare sRDA to sPLS, Both methods are applied on a
dataset with 150 patients and 200 mrna and 142 protein measurement.

5 components are extracted with enforcing 25 non-zero variables (no
optimization procedure). 


<a id="orgfd491ff"></a>

## Overlapping variables

    set.seed(100)
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
    
    #loadings.star and loadngs are the same for pls?
    res_spls$loadings$X[,1] == res_spls$loadings.star[[1]][,1]
    
    # Compare the variables selection
    res_sRDA <- get_nonzero_variables(res_object = res_sRDA)
    res_spls <- get_nonzero_variables(res_object = res_spls)
    
    #res_sRDA$nz_loading_names[["X"]]
    #res_spls$nz_loading_names[["X"]]
    
    plot_cus(1:ncomp,get_nr_of_common_components(res_sRDA, res_spls, 
    				       ncomp = ncomp), 
       ylab = "Nr of common variables", xlab = "Component", 
       ylim = c(0,nr_nonz))

<./plots/fig_overlap.pdf>

The variables selected in the first component
are identical for PLS and RDA, there are 5 overlapping variables in
the second component, and there are no overlapping variables in
consecutive latent variables (Figure [fig_overlap_sRDA_sPLS](#fig_overlap_sRDA_sPLS)).


<a id="orgea94e63"></a>

## Benchmark correlation

    # we can look at explained variances, they are about the 
    explained_Y <- cbind(res_spls$explained_variance$Y,
    res_sRDA$explained_variance$Y)
    
    colnames(explained_Y) <- c("Explined by sPLS",
    "Explained by sRDA")
    
    explained_X <- cbind(res_spls$explained_variance$X,
    res_sRDA$explained_variance$X)
    
    colnames(explained_X) <- c("Explined by sPLS",
    "Explained by sRDA")
    
    total_in_Y <- apply(explained_Y,2,sum)
    explained_Y <- rbind(explained_Y, total_in_Y)
    
    total_in_X <- apply(explained_X,2,sum)
    explained_X <- rbind(explained_X, total_in_X)
    
    round(explained_Y,2)

<table border="2" cellspacing="0" cellpadding="6" rules="groups" frame="hsides">


<colgroup>
<col  class="org-left" />

<col  class="org-right" />

<col  class="org-right" />
</colgroup>
<thead>
<tr>
<th scope="col" class="org-left">&#xa0;</th>
<th scope="col" class="org-right">Explined by sPLS</th>
<th scope="col" class="org-right">Explained by sRDA</th>
</tr>
</thead>

<tbody>
<tr>
<td class="org-left">comp 1</td>
<td class="org-right">0.13</td>
<td class="org-right">0.18</td>
</tr>


<tr>
<td class="org-left">comp 2</td>
<td class="org-right">0.12</td>
<td class="org-right">0.04</td>
</tr>


<tr>
<td class="org-left">comp 3</td>
<td class="org-right">0.06</td>
<td class="org-right">0.02</td>
</tr>


<tr>
<td class="org-left">comp 4</td>
<td class="org-right">0.07</td>
<td class="org-right">0.04</td>
</tr>


<tr>
<td class="org-left">comp 5</td>
<td class="org-right">0.04</td>
<td class="org-right">0.02</td>
</tr>


<tr>
<td class="org-left">total<sub>in</sub><sub>Y</sub></td>
<td class="org-right">0.41</td>
<td class="org-right">0.29</td>
</tr>
</tbody>
</table>

    
    round(explained_X,2)

<table border="2" cellspacing="0" cellpadding="6" rules="groups" frame="hsides">


<colgroup>
<col  class="org-left" />

<col  class="org-right" />

<col  class="org-right" />
</colgroup>
<thead>
<tr>
<th scope="col" class="org-left">&#xa0;</th>
<th scope="col" class="org-right">Explined by sPLS</th>
<th scope="col" class="org-right">Explained by sRDA</th>
</tr>
</thead>

<tbody>
<tr>
<td class="org-left">comp 1</td>
<td class="org-right">0.17</td>
<td class="org-right">0.17</td>
</tr>


<tr>
<td class="org-left">comp 2</td>
<td class="org-right">0.12</td>
<td class="org-right">0.12</td>
</tr>


<tr>
<td class="org-left">comp 3</td>
<td class="org-right">0.04</td>
<td class="org-right">0.05</td>
</tr>


<tr>
<td class="org-left">comp 4</td>
<td class="org-right">0.05</td>
<td class="org-right">0.04</td>
</tr>


<tr>
<td class="org-left">comp 5</td>
<td class="org-right">0.05</td>
<td class="org-right">0.03</td>
</tr>


<tr>
<td class="org-left">total<sub>in</sub><sub>X</sub></td>
<td class="org-right">0.43</td>
<td class="org-right">0.41</td>
</tr>
</tbody>
</table>


<a id="org6fd3c3b"></a>

# References


# Bibliography
<a id="legendre12_numer"></a>[legendre12_numer] Legendre, Numerical ecology, Elsevier (2012). [↩](#59e1f790b3b37a579e8a75fcebe25b6a)

<a id="Rohart_2017"></a>[Rohart_2017] Rohart, Gautier, Singh, & Lê Cao, mixOmics: an R package for ‘omics feature selection  and multiple data integration, <i></i>,  (2017). <a href="http://dx.doi.org/10.1101/108597">link</a>. <a href="http://dx.doi.org/10.1101/108597">doi</a>. [↩](#96148de61680d5941b2957b7ed478649)

<a id="wollenberg77_redun_analy_alter_canon_correl_analy"></a>[wollenberg77_redun_analy_alter_canon_correl_analy] Arnold van den Wollenberg, Redundancy Analysis an Alternative for Canonical  Correlation Analysis, <i>Psychometrika</i>, <b>42(2)</b>, 207-219 (1977). <a href="https://doi.org/10.1007/bf02294050">link</a>. <a href="http://dx.doi.org/10.1007/bf02294050">doi</a>. [↩](#789070d7c6e50336d15a21f274c89c09)

<a id="csala17_spars_redun_analy_high_dimen"></a>[csala17_spars_redun_analy_high_dimen] Attila Csala, Frans P J M Voorbraak, Aeilko H, Zwinderman & Michel H Hof, Sparse Redundancy Analysis of High-Dimensional  Genetic and Genomic Data, <i>Bioinformatics</i>, <b>33(20)</b>, 3228-3234 (2017). <a href="https://doi.org/10.1093/bioinformatics/btx374">link</a>. <a href="http://dx.doi.org/10.1093/bioinformatics/btx374">doi</a>. [↩](#294ead7176458d8a872380daa50f5068)

