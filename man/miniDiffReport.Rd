\name{miniDiffReport}
\alias{miniDiffReport}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Get XCMS-style diffReport from X13CMS data.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Generates conventional XCMS-like diff report comparing peaks between two sample classes, e.g. "control", and "treatment" from an xcmsSet object configured for X13CMS (i.e. one in which the sample classes have to be demarcated as unlabeled and labeled)
}
\usage{
miniDiffReport(xcmsSet, class1sampNames, class2sampNames, varEq = FALSE, intChoice)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xcmsSet}{
%%     ~~Describe \code{xcmsSet} here~~
XCMS object
}
  \item{class1sampNames}{
%%     ~~Describe \code{class1sampNames} here~~
The subset of rownames(xcmsSet@phenoData) corresponding to names of the unlabeled samples of class 1 (e.g. "control").
}
  \item{class2sampNames}{
%%     ~~Describe \code{class2sampNames} here~~
The subset of rownames(xcmsSet@phenoData) corresponding to names of the unlabeled samples of class 2 (e.g. "perturb").
}
  \item{varEq}{
%%     ~~Describe \code{varEq} here~~
Boolean indicating whether to assume that absolute ion intensities of each peak in both sample classes are drawn from distributions with equal variance. Defaults to FALSE.
}
  \item{intChoice}{
%%     ~~Describe \code{intChoice} here~~
one of "maxo", "into", or "intb"--the choice of which peak intensity measurement to use from the XCMS object.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
A dataframe whose rows correspond to peak groups identified in the xcmsSet object and whose columns are:
- pvalue: p-value of Welch's t-test comparing absolute intensities of the peak in the two conditions
- foldChange: change in mean peak intensity in condition 2 compared to conditions 1
- mzmed: median m/z of peak
- rtmed: median retention time of peak
- means1: mean peak intensity in condition 1
- means2: mean peak intensity in condition 2
- peakIntensities1: all peak intensities in the group in condition 1
- peakIntensities2: all peak intensities in the group in condition 2
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
Xiaojing Huang
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
## Not run:

## diffReport = miniDiffReport(xcmsSet, class1sampNames, class2sampNames, varEq = FALSE, intChoice = "intb") 

## From the example given in getIsoLabelReport(), the variables class1sampNames and class2sampNames would be c("control_unlabeled1", "control_unlabeled2"), and c("perturb_unlabeled1", "perturb_unlabeled2"), respectively.
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
