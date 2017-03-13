\name{plotLabelReport}
\alias{plotLabelReport}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Plot isotopologue groups from labeling report
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Generates bar plots of ion intensity distribution among the peaks of an isotopologue group
}
\usage{
plotLabelReport(isoLabelReport, intOption = "rel", classes, labeledSamples, outputfile)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{isoLabelReport}{
%%     ~~Describe \code{isoLabelReport} here~~
Output of getIsoLabelReport().
}
  \item{intOption}{
%%     ~~Describe \code{intOption} here~~
Choice of plotting relative ("rel") or absolute ("abs") ion intensities. Defaults to "rel".
}
  \item{classes}{
Character vector designating whether each sample is unlabeled (e.g. "C12") or labeled (e.g. "C13")
}
  \item{labeledSamples}{
%%     ~~Describe \code{labeledSamples} here~~
character variable designating labeled samples (e.g. "C13")
}
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
Pdf file containing 12 isotopologue group plots per page.
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
## plotLabelReport(labelsControl, "rel", "labelsControlPlots.pdf")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
