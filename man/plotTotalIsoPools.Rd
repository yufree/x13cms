\name{plotTotalIsoPools}
\alias{plotTotalIsoPools}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Plot output of getIsoDiffReport()
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Generates bar plots comparing isotopologue distributions in labeled samples of two different sample classes. Use absolute ion intensities to show distributions and compares total pool sizes of isotopologue group between sample classes.
}
\usage{
plotTotalIsoPools(isoDiffReport, xcmsSet, intChoice, sampleNames1, sampleNames2, labelReport1, labelReport2, classes1, classes2, labeledSamples, outputfile)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{isoDiffReport}{
%%     ~~Describe \code{isoDiffReport} here~~
Output of getIsoDiffReport().
}
  \item{xcmsSet}{xcmsSet object containing grouped and retention-time-aligned peaks (i.e. after calling group() and retcor() in XCMS)}
  \item{intChoice}{"intChoice" from getIsoLabelReport()}
  \item{sampleNames1}{"sampleNames" from getIsoLabelReport() called for condition 1}
  \item{sampleNames2}{"sampleNames" from getIsoLabelReport() called for condition 2}
  \item{labelReport1}{output of getIsoLabelReport() for condition 1}
  \item{labelReport2}{output of getIsoLabelReport() for condition 2}
  \item{classes1}{character vector designating whether each sample in labelsReport1 is unlabeled or labeled}
  \item{classes2}{character vector designating whether each sample in labelsReport2 is unlabeled or labeled}
  \item{labeledSamples}{character variable designating labeled samples (e.g. "C13")}
  \item{outputfile}{
%%     ~~Describe \code{outputfile} here~~
Name of pdf file to which plots are drawn.
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
%% ...
Pdf file containing 12 comparison plots per page.
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
## plotTotalPools(isoDiffReport, "totalPools.pdf") 
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
