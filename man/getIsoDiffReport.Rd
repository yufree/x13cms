\name{getIsoDiffReport}
\alias{getIsoDiffReport}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Compare isotope labeling reports 
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This routine matches up label-enriched isotopologue groups from each biological condition and determines whether the labeling patterns in each group are different between the conditions. To print the list to a tab-delimited file, use \code{\link{printIsoListOutputs}}.
}
\usage{
getIsoDiffReport(labelsData1, labelsData2, condition1, condition2, classes1, classes2, labeledSamples, varEq, singleSample)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{labelsData1}{
%%     ~~Describe \code{labelsData1} here~~
isoLabelReport for samples from conditions 1
}
  \item{labelsData2}{
%%     ~~Describe \code{labelsData2} here~~
isoLabelReport for samples from conditions 2
}
  \item{condition1}{
%%     ~~Describe \code{labelsData1} here~~
character variable labeling condition 1, e.g. "control"
}
  \item{condition2}{
%%     ~~Describe \code{labelsData2} here~~
character variable labeling condition 2, e.g. "perturb"
}
  \item{classes1}{
%%     ~~Describe \code{classes1} here~~
character vector designating whether each sample in labelsData1 is unlabeled or labeled
}
  \item{classes2}{
%%     ~~Describe \code{classes2} here~~
character vector designating whether each sample in labelsData2 is unlabeled or labeled
}
  \item{labeledSamples}{
%%     ~~Describe \code{labeledSamples} here~~
character variable designating labeled samples (e.g. "C13")
}
  \item{varEq}{
%%     ~~Describe \code{varEq} here~~
Boolean indicating whether to assume that relative isotopologue intensities in each condition are drawn from distributions with equal variance. Defaults to FALSE.
}
  \item{singleSample}{
%%     ~~Describe \code{varEq} here~~
Boolean indicating whether only single samples were used to generate the labeling reports. Defaults to FALSE
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1}{Description of 'comp1'}
%%  \item{comp2}{Description of 'comp2'}
%% ...
An isoDiffReport consisting of the following lists describing isotopologue groups:
\item{condition1}{m/z of all isotopologues in a group in samples from condition 1; empty if group is found to be enriched for label only in condition 2}
\item{condition2}{m/z of all isotopologues in a group in samples from condition 2; empty if group is found to be enriched for label only in condition 1}
\item{rt1}{Retention time of each isotopologue in condition 1}
\item{rt2}{Retention time of each isotopologue in condition 2}
\item{relInts1U}{Mean relative intensities of each isotopologue in the samples of condition 1 treated with unlabeled precursor}
\item{relInts1L}{Mean relative intensities of each isotopologue in the samples of condition 1 treated with labeled precursor}
\item{relInts2U}{Mean relative intensities of each isotopologue in the samples of condition 2 treated with unlabeled precursor}
\item{relInts2L}{Mean relative intensities of each isotopologue in the samples of condition 2 treated with labeled precursor}
\item{p_value}{p-values from Welch's t-tests comparing relative intensities of each isotopologue in labeled samples of condition 1 vs. those of condition 2}
\item{sdRelInts1L}{std dev of relative intensity of each isotopologue peak in labeled samples of condition 1}
\item{sdRelInts2L}{std dev of relative intensity of each isotopologue peak in labeled samples of condition 2}
\item{absInts1L}{mean absolute ion intensity of each isotopologue peak in labeled samples of condition 1}
\item{absInts2L}{mean absolute ion intensity of each isotopologue peak in labeled samples of condition 2}
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

%% ~Make other sections like Warning with \section{Warning}{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
## Not run: 
## isoDiffReport = getIsoDiffReport(labelsControl, labelsPerturb, condition1 = "control", condition2 = "perturb", classesControl = c("C12", "C12", "C13", "C13"), classesPerturb = c("C12", "C12", "C13", "C13"), labeledSamples = "C13", varEq = FALSE)
##
## This command generates the isoDiffReport for the example given for \code{\link{getIsoLabelReport}}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1}
\keyword{ ~kwd2}% __ONLY ONE__ keyword per line
