\name{getIsoLabelReport}
\alias{getIsoLabelReport}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Generate isotope labeling report
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function searches the groups table of the input xcmsSet object for peaks that are potential isotopologues of one another. It determines whether each isotopologue group identified is significantly enriched for the added isotope label. The output is a list of all enriched groups, along with information about each one's labeling pattern. To print the list to a tab-delimited file, use printIsoListOutputs.
}
\usage{
getIsoLabelReport(xcmsSet, sampleNames, unlabeledSamples, labeledSamples, isotopeMassDiff, RTwindow, ppm, massOfLabeledAtom, noiseCutoff, intChoice = "intb", varEq = FALSE, alpha, singleSample = FALSE, compareOnlyDistros = FALSE, monotonicityTol = FALSE, enrichTol = 0.1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xcmsSet}{
%%     ~~Describe \code{xcmsSet} here~~
xcmsSet object containing grouped and retention-time-aligned peaks (i.e. after calling group() and retcor() in XCMS); sample classes in the @phenoData slot should be designated with the arguments "unlabeledSamples" and "labeledSamples"
}
  \item{sampleNames}{
%%     ~~Describe \code{sampleNames} here~~
character vector of names of the unlabeled and labeled samples for which the labeling report should be generated; taken from rownames of xcmsSet@phenoData
}
  \item{unlabeledSamples}{
%%     ~~Describe \code{unlabeledSamples} here~~
character variable designating unlabeled samples (e.g. "C12"). Must match corresponding entries of @phenoData.
}
  \item{labeledSamples}{
%%     ~~Describe \code{labeledSamples} here~~
character variable designating labeled samples (e.g. "C13"). Must match corresponding entries of @phenoData.
}
  \item{isotopeMassDiff}{
%%     ~~Describe \code{isotopeMassDiff} here~~
difference in mass between labeled and unlabeled atom (e.g. 1.00335 for C13) 
}
  \item{RTwindow}{
%%     ~~Describe \code{RTwindow} here~~
retention time window in which all peaks are considered to be co-eluting, in seconds
}
  \item{ppm}{
%%     ~~Describe \code{ppm} here~~
ppm allowance for deviation of peaks within an isotopologue group from expected m/z; in practice this should be set higher than ppm tol used for peak-picking (e.g. 20 for a 5 ppm instrument) to ensure that all isotopologues are captured
}
  \item{massOfLabeledAtom}{
%%     ~~Describe \code{massOfLabeledAtom} here~~
e.g. 12.0000 for C12
}
  \item{noiseCutoff}{
%%     ~~Describe \code{noiseCutoff} here~~
ion intensity cutoff below which a peak is considered noise
}
  \item{intChoice}{
%%     ~~Describe \code{intChoice} here~~
one of "maxo", "into", or "intb"--the choice of which peak intensity measurement to use from the XCMS object. Defaults to "intb".
}
  \item{varEq}{
%%     ~~Describe \code{varEq} here~~
Boolean indicating whether to assume equal variance on the peak intensities in the unlabeled and labeled samples. Defaults to FALSE.
}
  \item{alpha}{
%%     ~~Describe \code{alpha} here~~
p-value cutoff for calling significance of label enrichment
}
  \item{singleSample}{
%%     ~~Describe \code{asingleSample} here~~
Boolean indicating whether only single replicates exist for unlabeled and labeled samples. Defaults to FALSE.
}
  \item{compareOnlyDistros}{
%%     ~~Describe \code{asingleSample} here~~
Boolean indicating whether to look at distributions with intent of calling enrichment over unlabeled (FALSE, default) or to compare distributions of label between two labeled samples of different classes (TRUE).
}
  \item{monotonicityTol}{
Tolerance parameter used to enforce expected ion intensity pattern (i.e. monotonic decrease from M0 to Mn) in unlabeled samples; a low value closer to 0 enforces stricter monotonicity; default is to not enforce monotonicity (monotonicityTol = FALSE) due to potential carryover between samples}
  \item{enrichTol}{
tolerance parameter for enforcing enrichment of higher isotopologues in labeled samples; a value of 0 enforces strict requirement for enrichment of higher isotopologues to be higher in labeled samples
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
An isoLabelReport consisting of the following lists:
\item{compound}{m/z of the base (unlabeled) isotopologue; value to be used for searching metabolite databases for putative identities of the group}
\item{isotologue}{List of m/z of all isotopologues in the group, including base}
\item{groupID}{ID number of isotopologues in the xcmsSet object's @groups slot; used for plotting EICs of individual isotopologues}
\item{rt}{Retention time of each isotopologue}
\item{meanAbsIntU}{Mean absolute intensity of each isotopologue in the samples treated with unlabeled precursor}
\item{totalAbsU}{Total ion intensity for all isotopologues in a group}
\item{cvTotalU}{coefficient of variation of total ion intensity of isotopologue group in unlabeled samples}
\item{meanAbsIntL}{Mean absolute intensity of each isotopologue in the samples treated with labeled precursor}
\item{totalAbsL}{Total ion intensity for all isotopologues in a group}
\item{cvTotalL}{coefficient of variation of total ion intensity of isotopologue group in unlabeled samples}
\item{meanRelIntU}{Mean relative intensity (fraction of total absolute intensity of entire isotopologue group) of each isotopologue in the samples treated with unlabeled precursor}
\item{meanRelIntL}{Mean relative intensity of each isotopologue in the samples treated with labeled precursor}
\item{p_value}{p-values from Welch's t-tests comparing relative intensities of each isotopologue in unlabeled vs. labeled samples}
\item{enrichmentLvsU}{Ratio meanRelIntL: meanRelIntU  for each isotopologue}
\item{sdRelU}{std dev of relative intensity of each isotopologue peak in unlabeled samples}
\item{sdRelL}{std dev of relative intensity of each isotopologue peak in labeled samples}
\item{sampleData}{Absolute intensities of each isotopologue in every sample}

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
## labelsControl = getIsoLabelReport(xcmsSet = xs, sampleNames = c("control_unlabeled1", "control_unlabeled2", "control_labeled1", "control_labeled2"), unlabeledSamples = "C12", labeledSamples = "C13", isotopeMassDiff = 1.00335, RTwindow = 10, ppm = 20, massOfLabeledAtom = 12, noiseCutoff = 10000, intChoice = "intb", varEq = FALSE, alpha = 0.05, singleSample = FALSE, compareOnlyDistros = FALSE, monotonicityTol = 0.1)
## 
## labelsPerturb = getIsoLabelReport(xcmsSet = xs, sampleNames = c("perturb_unlabeled1", "perturb_unlabeled2", "perturb_labeled1", "perturb_labeled2"), unlabeledSamples = "C12", labeledSamples = "C13", isotopeMassDiff = 1.00335, RTwindow = 10, ppm = 20, massOfLabeledAtom = 12, noiseCutoff = 10000, intChoice = "intb", varEq = FALSE, alpha = 0.05, singleSample = FALSE, compareOnlyDistros = FALSE, monotonicityTol = 0.1)
##
## This experiment consists of 2 replicates of unlabeled samples and 2 of labeled samples in each biological condition, control or perturbed. The @phenoData slot of the xcmsSet object "xs" would have rownames = c("control_unlabeled1", "control_unlabeled2", "perturb_unlabeled1", "perturb_unlabeled2", "control_labeled1", "control_labeled2", "perturb_labeled1", "perturb_labeled2") and values = c("C12", "C12", "C12", "C12", "C13", "C13", "C13", "C13").
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1}
\keyword{ ~kwd2}% __ONLY ONE__ keyword per line
