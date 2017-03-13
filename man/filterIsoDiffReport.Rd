\name{filterIsoDiffReport}
\alias{filterIsoDiffReport}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Filter isoDiffReport 
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This routine filters an isoDiffReport to include only isotopologue groups that are differentially labeled between sample conditions and are not falsely unique due to suboptimal chromatography conditions. An additional, optional filter exists for the special case of 13C labeling experiments; it removes groups that consist only of M and M+1 peaks.
}
\usage{
filterIsoDiffReport(isoDiffReport, alpha, whichPeak = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{isoDiffReport}{
%%     ~~Describe \code{isoDiffReport} here~~
report to be filtered
}
  \item{alpha}{
%%     ~~Describe \code{alpha} here~~
p-value cutoff to call labeling pattern of isotoplogue group as different between sample conditions.
}
  \item{whichPeak}{
%%     ~~Describe \code{whichPeak} here~~
option for which isotopologue peaks within the group to consider when calling significance for differential labeling. Default = 1 indicates consideration of the relative intensities of the base peak only. '2' indicates consideration of all peaks and a call of significance if at least one of their relative intensities is significantly different between the two sample conditions.
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
isoDiffReport
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
## filterIsoDiffReport(isoDiffReport, alpha = 0.05, whichPeak = 1, is13C = FALSE)   
 
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
