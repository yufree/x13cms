require(xcms)
require(stats)
getPeakIntensities <- function(xcmsSet, sampleNames, 
    intChoice) {
    peakIntensities = groupval(xcmsSet, method = "medret", 
        value = intChoice)
    peakIntensities = peakIntensities[, match(sampleNames, 
        colnames(peakIntensities))]
    if (intChoice == "intb") {
        peakIntensities[is.na(peakIntensities)] = 0
    }
    return(peakIntensities)
}


#' %% ~~function to do ... ~~ Generate isotope labeling report
#' 
#' %% ~~ A concise (1-5 lines) description of what the function does. ~~ This
#' function searches the groups table of the input xcmsSet object for peaks
#' that are potential isotopologues of one another. It determines whether each
#' isotopologue group identified is significantly enriched for the added
#' isotope label. The output is a list of all enriched groups, along with
#' information about each one's labeling pattern. To print the list to a
#' tab-delimited file, use printIsoListOutputs.
#' 
#' %% ~~ If necessary, more details than the description above ~~
#' 
#' @param xcmsSet %% ~~Describe \code{xcmsSet} here~~ xcmsSet object containing
#' grouped and retention-time-aligned peaks (i.e. after calling group() and
#' retcor() in XCMS); sample classes in the @phenoData slot should be
#' designated with the arguments 'unlabeledSamples' and 'labeledSamples'
#' @param sampleNames %% ~~Describe \code{sampleNames} here~~ character vector
#' of names of the unlabeled and labeled samples for which the labeling report
#' should be generated; taken from rownames of xcmsSet@phenoData
#' @param unlabeledSamples %% ~~Describe \code{unlabeledSamples} here~~
#' character variable designating unlabeled samples (e.g. 'C12'). Must match
#' corresponding entries of @phenoData.
#' @param labeledSamples %% ~~Describe \code{labeledSamples} here~~ character
#' variable designating labeled samples (e.g. 'C13'). Must match corresponding
#' entries of @phenoData.
#' @param isotopeMassDiff %% ~~Describe \code{isotopeMassDiff} here~~
#' difference in mass between labeled and unlabeled atom (e.g. 1.00335 for C13)
#' @param RTwindow %% ~~Describe \code{RTwindow} here~~ retention time window
#' in which all peaks are considered to be co-eluting, in seconds
#' @param ppm %% ~~Describe \code{ppm} here~~ ppm allowance for deviation of
#' peaks within an isotopologue group from expected m/z; in practice this
#' should be set higher than ppm tol used for peak-picking (e.g. 20 for a 5 ppm
#' instrument) to ensure that all isotopologues are captured
#' @param massOfLabeledAtom %% ~~Describe \code{massOfLabeledAtom} here~~ e.g.
#' 12.0000 for C12
#' @param noiseCutoff %% ~~Describe \code{noiseCutoff} here~~ ion intensity
#' cutoff below which a peak is considered noise
#' @param intChoice %% ~~Describe \code{intChoice} here~~ one of 'maxo',
#' 'into', or 'intb'--the choice of which peak intensity measurement to use
#' from the XCMS object. Defaults to 'intb'.
#' @param varEq %% ~~Describe \code{varEq} here~~ Boolean indicating whether to
#' assume equal variance on the peak intensities in the unlabeled and labeled
#' samples. Defaults to FALSE.
#' @param alpha %% ~~Describe \code{alpha} here~~ p-value cutoff for calling
#' significance of label enrichment
#' @param singleSample %% ~~Describe \code{asingleSample} here~~ Boolean
#' indicating whether only single replicates exist for unlabeled and labeled
#' samples. Defaults to FALSE.
#' @param compareOnlyDistros %% ~~Describe \code{asingleSample} here~~ Boolean
#' indicating whether to look at distributions with intent of calling
#' enrichment over unlabeled (FALSE, default) or to compare distributions of
#' label between two labeled samples of different classes (TRUE).
#' @param monotonicityTol Tolerance parameter used to enforce expected ion
#' intensity pattern (i.e. monotonic decrease from M0 to Mn) in unlabeled
#' samples; a low value closer to 0 enforces stricter monotonicity; default is
#' to not enforce monotonicity (monotonicityTol = FALSE) due to potential
#' carryover between samples
#' @param enrichTol tolerance parameter for enforcing enrichment of higher
#' isotopologues in labeled samples; a value of 0 enforcesct requirement
#' for enrichment of higher isotopologues to be higher in labeled samples
#' @return %% ~Describe the value returned %% If it is a LIST, use %%
#' \item{comp1}{Description of 'comp1'} %% \item{comp2}{Description of 'comp2'}
#' %% ... An isoLabelReport consisting of the following lists:
#' \item{compound}{m/z of the base (unlabeled) isotopologue; value to be used
#' for searching metabolite databases for putative identities of the group}
#' \item{isotologue}{List of m/z of all isotopologues in the group, including
#' base} \item{groupID}{ID number of isotopologues in the xcmsSet object's
#' @groups slot; used for plotting EICs of individual isotopologues}
#' \item{rt}{Retention time of each isotopologue} \item{meanAbsIntU}{Mean
#' absolute intensity of each isotopologue in the samples treated with
#' unlabeled precursor} \item{totalAbsU}{Total ion intensity for all
#' isotopologues in a group} \item{cvTotalU}{coefficient of variation of total
#' ion intensity of isotopologue group in unlabeled samples}
#' \item{meanAbsIntL}{Mean absolute intensity of each isotopologue in the
#' samples treated with labeled precursor} \item{totalAbsL}{Total ion intensity
#' for all isotopologues in a group} \item{cvTotalL}{coefficient of variation
#' of total ion intensity of isotopologue group in unlabeled samples}
#' \item{meanRelIntU}{Mean relative intensity (fraction of total absolute
#' intensity of entire isotopologue group) of each isotopologue in the samples
#' treated with unlabeled precursor} \item{meanRelIntL}{Mean relative intensity
#' of each isotopologue in the samples treated with labeled precursor}
#' \item{p_value}{p-values from Welch's t-tests comparing relative intensities
#' of each isotopologue in unlabeled vs. labeled samples}
#' \item{enrichmentLvsU}{Ratio meanRelIntL: meanRelIntU for each isotopologue}
#' \item{sdRelU}{std dev of relative intensity of each isotopologue peak in
#' unlabeled samples} \item{sdRelL}{std dev of relative intensity of each
#' isotopologue peak in labeled samples} \item{sampleData}{Absolute intensities
#' of each isotopologue in every sample}
#' @note %% ~~further notes~~
#' @author %% ~~who you are~~ Xiaojing Huang
#' @seealso %% ~~objects to See Also as \code{\link{help}}, ~~~
#' @references %% ~put references to the literature/web site here ~
#' @keywords ~kwd1 ~kwd2
#' @examples
#' 
#' ## Not run: 
#' ## labelsControl = getIsoLabelReport(xcmsSet = xs, sampleNames = c('control_unlabeled1', 'control_unlabeled2', 'control_labeled1', 'control_labeled2'), unlabeledSamples = 'C12', labeledSamples = 'C13', isotopeMassDiff = 1.00335, RTwindow = 10, ppm = 20, massOfLabeledAtom = 12, noiseCutoff = 10000, intChoice = 'intb', varEq = FALSE, alpha = 0.05, singleSample = FALSE, compareOnlyDistros = FALSE, monotonicityTol = 0.1)
#' ## 
#' ## labelsPerturb = getIsoLabelReport(xcmsSet = xs, sampleNames = c('perturb_unlabeled1', 'perturb_unlabeled2', 'perturb_labeled1', 'perturb_labeled2'), unlabeledSamples = 'C12', labeledSamples = 'C13', isotopeMassDiff = 1.00335, RTwindow = 10, ppm = 20, massOfLabeledAtom = 12, noiseCutoff = 10000, intChoice = 'intb', varEq = FALSE, alpha = 0.05, singleSample = FALSE, compareOnlyDistros = FALSE, monotonicityTol = 0.1)
#' ##
#' ## This experiment consists of 2 replicates of unlabeled samples and 2 of labeled samples in each biological condition, control or perturbed. The @phenoData slot of the xcmsSet object 'xs' would have rownames = c('control_unlabeled1', 'control_unlabeled2', 'perturb_unlabeled1', 'perturb_unlabeled2', 'control_labeled1', 'control_labeled2', 'perturb_labeled1', 'perturb_labeled2') and values = c('C12', 'C12', 'C12', 'C12', 'C13', 'C13', 'C13', 'C13').
#' 
getIsoLabelReport <- function(xcmsSet, sampleNames, 
    unlabeledSamples, labeledSamples, isotopeMassDiff, 
    RTwindow, ppm, massOfLabeledAtom, noiseCutoff, 
    intChoice = "intb", varEq = FALSE, alpha, singleSample = FALSE, 
    compareOnlyDistros = FALSE, monotonicityTol = FALSE, 
    enrichTol = 0.1) {
    groups = data.frame(xcmsSet@groups)
    peakIntensities = getPeakIntensities(xcmsSet, 
        sampleNames, intChoice)
    peakIntensities = peakIntensities[order(groups$mzmed), 
        ]
    groups = groups[order(groups$mzmed), ]
    groupRTs = groups$rtmed
    groupMzs = groups$mzmed
    groupIDs = as.numeric(rownames(groups))
    nGroups = length(groupMzs)
    classes = as.character(xcmsSet@phenoData[match(sampleNames, 
        rownames(xcmsSet@phenoData)), ])
    numSamples = length(classes)
    intensities1 = peakIntensities[, which(classes == 
        unlabeledSamples), drop = FALSE]
    intensities2 = peakIntensities[, which(classes == 
        labeledSamples), drop = FALSE]
    iMD = isotopeMassDiff
    base = list()
    labeled = list()
    basePeak = list()
    labeledPeak = list()
    groupIndicesByRT = order(groupRTs)
    orderedGroupRTs = groupRTs[groupIndicesByRT]
    for (i in 1:nGroups) {
        binI = groupIndicesByRT[orderedGroupRTs - 
            orderedGroupRTs[i] >= 0 & orderedGroupRTs - 
            orderedGroupRTs[i] <= RTwindow]
        bin = groups[binI, ]
        binSize = length(binI)
        I = groupIndicesByRT[i]
        if (binSize > 0) {
            for (j in 1:binSize) {
                if (groups$mzmed[I] < bin$mzmed[j]) {
                  a = I
                  b = binI[j]
                } else {
                  a = binI[j]
                  b = I
                }
                delta = (groupMzs[b] - groupMzs[a])/iMD
                DELTA = round(delta)
                if (DELTA == 0) {
                  next
                }
                if (delta <= DELTA * (1 + ppm/1e+06) + 
                  (groupMzs[a] * ppm/1e+06)/(iMD * 
                    (1 - ppm/1e+06)) & delta >= DELTA * 
                  (1 - ppm/1e+06) - (groupMzs[a] * 
                  ppm/1e+06)/(iMD * (1 + ppm/1e+06))) {
                  if (DELTA * massOfLabeledAtom >= 
                    groupMzs[a]) {
                    next
                  }
                  if (mean(intensities1[b, ]) > mean(intensities1[a, 
                    ]) & !compareOnlyDistros) {
                    next
                  }
                  if (all(intensities1[a, ] == 0) & 
                    all(intensities2[a, ] == 0)) {
                    next
                  }
                  if (all(intensities1[b, ] == 0) & 
                    all(intensities2[b, ] == 0)) {
                    next
                  }
                  base = c(base, a)
                  labeled = c(labeled, b)
                  basePeak = c(basePeak, groupMzs[a])
                  labeledPeak = c(labeledPeak, groupMzs[b])
                }
            }
        }
    }
    labelsMatrix = as.matrix(cbind(unlist(base), unlist(labeled), 
        unlist(basePeak), unlist(labeledPeak)))
    labelsMatrix = labelsMatrix[order(labelsMatrix[, 
        3], labelsMatrix[, 4]), ]
    numPutativeLabels = dim(labelsMatrix)[1]
    basePeaks = unique(labelsMatrix[, 1])
    numLabeledPeaks = length(basePeaks)
    outtakes = list()
    for (i in 1:numPutativeLabels) {
        B = labelsMatrix[, 2] == labelsMatrix[i, 1]
        A = labelsMatrix[B, 1]
        C = which(labelsMatrix[, 1] %in% A)
        if (any(labelsMatrix[C, 2] == labelsMatrix[i, 
            2])) {
            outtakes = c(outtakes, i)
            next
        }
        if (i < numPutativeLabels) {
            A = (i + 1):numPutativeLabels
            idx = any(labelsMatrix[A, 1] == labelsMatrix[i, 
                1] & labelsMatrix[A, 2] == labelsMatrix[i, 
                2])
            if (idx) {
                outtakes = c(outtakes, i)
            }
        }
    }
    outtakes = unlist(outtakes)
    labelsMatrix = labelsMatrix[-outtakes, ]
    numPutativeLabels = dim(labelsMatrix)[1]
    basePeaks = unique(labelsMatrix[, 1])
    numLabeledPeaks = length(basePeaks)
    base = list()
    mz = list()
    ID = list()
    RT = list()
    absInt1 = list()
    absInt2 = list()
    relInt1 = list()
    relInt2 = list()
    totInt1 = list()
    totInt2 = list()
    CVabsInt1 = list()
    CVabsInt2 = list()
    SDrelInt1 = list()
    SDrelInt2 = list()
    foldEnrichment = list()
    pvalues = list()
    sampleIntensities = list()
    j = 1
    for (i in 1:numLabeledPeaks) {
        a = basePeaks[i]
        baseIntensities = c(intensities1[a, ], intensities2[a, 
            ])
        isotopologues = list()
        IDs = list()
        RTs = list()
        numisotopologues = 0
        k = j
        while (k <= numPutativeLabels){
            if(labelsMatrix[k, 1] != a){
            break
        }
            isotopologues = c(isotopologues, groupMzs[labelsMatrix[k, 
                2]])
            IDs = c(IDs, groupIDs[labelsMatrix[k, 
                2]])
            RTs = c(RTs, groupRTs[labelsMatrix[k, 
                2]])
            numisotopologues = numisotopologues + 
                1
            k = k + 1
        }
        isotopologues = unlist(isotopologues)
        IDs = unlist(IDs)
        RTs = unlist(RTs)
        if (mean(intensities1[a, ]) < noiseCutoff) {
            j = k
            next
        }
        abs1 = list()
        abs2 = list()
        labeledIntensities = matrix(rep(0, numisotopologues * 
            numSamples), nrow = numisotopologues, 
            ncol = numSamples)
        for (l in 1:numisotopologues) {
            b = labelsMatrix[j + l - 1, 2]
            labeledIntensities[l, ] = cbind(intensities1[b, 
                , drop = FALSE], intensities2[b, , 
                drop = FALSE])
            abs1 = c(abs1, mean(intensities1[b, ]))
            abs2 = c(abs2, mean(intensities2[b, ]))
        }
        abs1 = unlist(abs1)
        abs2 = unlist(abs2)
        if (numisotopologues != length(unique(round(isotopologues)))) {
            M0 = round(groupMzs[a])
            isos = round(isotopologues)
            reduced = unique(isos)
            numUniqIsos = length(reduced)
            outtakes = list()
            for (r in 1:numUniqIsos) {
                q = which(isos == reduced[r])
                if (length(q) > 1) {
                  massdefect = iMD * (reduced[r] - 
                    M0)
                  delta = abs(groupMzs[IDs[q]] - groupMzs[a] - 
                    massdefect)
                  outtakes = c(outtakes, q[which(delta != 
                    min(delta))])
                }
            }
            if (length(outtakes) > 0) {
                outtakes = unlist(outtakes)
                isotopologues = isotopologues[-outtakes]
                IDs = IDs[-outtakes]
                RTs = RTs[-outtakes]
                numisotopologues = length(isotopologues)
                abs1 = abs1[-outtakes]
                abs2 = abs2[-outtakes]
                labeledIntensities = labeledIntensities[-outtakes, 
                  , drop = FALSE]
            }
        }
        if (!compareOnlyDistros & monotonicityTol) {
            meanMprevUL = mean(intensities1[a, ])
            outtakes = list()
            for (l in 1:numisotopologues) {
                if (l == 1 & abs1[l] > (1 + monotonicityTol) * 
                  meanMprevUL) {
                  outtakes = c(outtakes, l)
                } else if (l > 1 & abs1[l] > (1 + monotonicityTol) * 
                  meanMprevUL & round(isotopologues[l] - 
                  isotopologues[l - 1]) > 1) {
                  outtakes = c(outtakes, l)
                } else {
                  meanMprevUL = abs1[l]
                }
            }
            outtakes = unlist(outtakes)
            if (length(outtakes) > 0) {
                abs1 = abs1[-outtakes]
                abs2 = abs2[-outtakes]
                labeledIntensities = labeledIntensities[-outtakes, 
                  , drop = FALSE]
                isotopologues = isotopologues[-outtakes]
                IDs = IDs[-outtakes]
                RTs = RTs[-outtakes]
            }
        }
        isotopologues = c(groupMzs[a], isotopologues)
        IDs = c(groupIDs[a], IDs)
        RTs = c(groupRTs[a], RTs)
        allIntensities = rbind(baseIntensities, labeledIntensities)
        abs1 = c(mean(intensities1[a, ]), abs1)
        abs2 = c(mean(intensities2[a, ]), abs2)
        numisotopologues = length(isotopologues)
        sumIntensities = colSums(allIntensities)
        tot1 = mean(sumIntensities[1:dim(intensities1)[2]])
        tot2 = mean(sumIntensities[(dim(intensities1)[2] + 
            1):numSamples])
        cv1 = sd(sumIntensities[1:dim(intensities1)[2]])/tot1
        cv2 = sd(sumIntensities[(dim(intensities1)[2] + 
            1):numSamples])/tot2
        groupIntensities = allIntensities/matrix(rep(sumIntensities, 
            numisotopologues), nrow = numisotopologues, 
            byrow = TRUE)
        gI1 = groupIntensities[, 1:dim(intensities1)[2], 
            drop = FALSE]
        gI2 = groupIntensities[, (dim(intensities1)[2] + 
            1):numSamples, drop = FALSE]
        gI1 = gI1[, colSums(is.na(gI1)) == 0, drop = FALSE]
        gI2 = gI2[, colSums(is.na(gI2)) == 0, drop = FALSE]
        if (dim(gI1)[2] < dim(intensities1)[2]/2 || 
            dim(gI2)[2] < dim(intensities2)[2]/2) {
            j = k
            next
        }
        rel1 = rowMeans(gI1)
        rel2 = rowMeans(gI2)
        sd1 = apply(gI1, 1, sd)
        sd2 = apply(gI2, 1, sd)
        enrichRatios = rel2/rel1
        if (!compareOnlyDistros) {
            if (enrichRatios[1] > (1 + enrichTol)) {
                j = k
                next
            }
        }
        if (!singleSample) {
            pvalue = list()
            for (l in 1:numisotopologues) {
                if (all(gI1[l, ] == 1) & all(gI2[l, 
                  ] == 0) || is.infinite(enrichRatios[l])) {
                  pvalue = c(pvalue, 0)
                } else {
                  T = try(t.test(gI1[l, ], gI2[l, 
                    ], var.equal = varEq), silent = TRUE)
                  if (class(T) == "try-error") {
                    pvalue = c(pvalue, 1)
                    break
                  } else {
                    pvalue = c(pvalue, T$p.value)
                  }
                }
            }
            if (any(unlist(pvalue) < alpha) & !any(unlist(pvalue) == 
                1)) {
                base = c(base, groupMzs[a])
                mz = c(mz, list(isotopologues))
                ID = c(ID, list(IDs))
                RT = c(RT, list(RTs))
                absInt1 = c(absInt1, list(abs1))
                absInt2 = c(absInt2, list(abs2))
                relInt1 = c(relInt1, list(rel1))
                relInt2 = c(relInt2, list(rel2))
                CVabsInt1 = c(CVabsInt1, cv1)
                CVabsInt2 = c(CVabsInt2, cv2)
                totInt1 = c(totInt1, tot1)
                totInt2 = c(totInt2, tot2)
                SDrelInt1 = c(SDrelInt1, list(sd1))
                SDrelInt2 = c(SDrelInt2, list(sd2))
                foldEnrichment = c(foldEnrichment, 
                  list(enrichRatios))
                pvalues = c(pvalues, list(unlist(pvalue)))
                sampleIntensities = c(sampleIntensities, 
                  list(allIntensities))
            }
        } else {
            deltaSpec = sum(abs(rel1[1:numisotopologues - 
                1] - rel2[1:numisotopologues - 1]))
            base = c(base, groupMzs[a])
            mz = c(mz, list(isotopologues))
            ID = c(ID, list(IDs))
            RT = c(RT, list(RTs))
            absInt1 = c(absInt1, list(abs1))
            absInt2 = c(absInt2, list(abs2))
            relInt1 = c(relInt1, list(rel1))
            relInt2 = c(relInt2, list(rel2))
            CVabsInt1 = c(CVabsInt1, cv1)
            CVabsInt2 = c(CVabsInt2, cv2)
            totInt1 = c(totInt1, tot1)
            totInt2 = c(totInt2, tot2)
            SDrelInt1 = c(SDrelInt1, list(sd1))
            SDrelInt2 = c(SDrelInt2, list(sd2))
            foldEnrichment = c(foldEnrichment, list(enrichRatios))
            pvalues = c(pvalues, deltaSpec)
            sampleIntensities = c(sampleIntensities, 
                list(allIntensities))
        }
        j = k
    }
    labelsData = list(compound = base, isotopologue = mz, 
        groupID = ID, rt = RT, meanAbsU = absInt1, 
        totalAbsU = totInt1, cvTotalU = CVabsInt1, 
        meanAbsL = absInt2, totalAbsL = totInt2, cvTotalL = CVabsInt2, 
        meanRelU = relInt1, meanRelL = relInt2, p_value = pvalues, 
        enrichmentLvsU = foldEnrichment, sdRelU = SDrelInt1, 
        sdRelL = SDrelInt2, sampleData = sampleIntensities)
    return(labelsData)
}


#' %% ~~function to do ... ~~ Compare isotope labeling reports
#' 
#' %% ~~ A concise (1-5 lines) description of what the function does. ~~ This
#' routine matches up label-enriched isotopologue groups from each biological
#' condition and determines whether the labeling patterns in each group are
#' different between the conditions. To print the list to a tab-delimited file,
#' use \code{\link{printIsoListOutputs}}.
#' 
#' %% ~~ If necessary, more details than the description above ~~
#' 
#' @param labelsData1 %% ~~Describe \code{labelsData1} here~~ isoLabelReport
#' for samples from conditions 1
#' @param labelsData2 %% ~~Describe \code{labelsData2} here~~ isoLabelReport
#' for samples from conditions 2
#' @param condition1 %% ~~Describe \code{labelsData1} here~~ character variable
#' labeling condition 1, e.g. 'control'
#' @param condition2 %% ~~Describe \code{labelsData2} here~~ character variable
#' labeling condition 2, e.g. 'perturb'
#' @param classes1 %% ~~Describe \code{classes1} here~~ character vector
#' designating whether each sample in labelsData1 is unlabeled or labeled
#' @param classes2 %% ~~Describe \code{classes2} here~~ character vector
#' designating whether each sample in labelsData2 is unlabeled or labeled
#' @param labeledSamples %% ~~Describe \code{labeledSamples} here~~ character
#' variable designating labeled samples (e.g. 'C13')
#' @param varEq %% ~~Describe \code{varEq} here~~ Boolean indicating whether to
#' assume that relative isotopologue intensities in each condition are drawn
#' from distributions with equal variance. Defaults to FALSE.
#' @param singleSample %% ~~Describe \code{varEq} here~~ Boolean indicating
#' whether only single samples were used to generate the labeling reports.
#' Defaults to FALSE
#' @return %% ~Describe the value returned %% If it is a LIST, use %%
#' \item{comp1}{Description of 'comp1'} %% \item{comp2}{Description of 'comp2'}
#' %% ... An isoDiffReport consisting of the following lists describing
#' isotopologue groups: \item{condition1}{m/z of all isotopologues in a group
#' in samples from condition 1; empty if group is found to be enriched for
#' label only in condition 2} \item{condition2}{m/z of all isotopologues in a
#' group in samples from condition 2; empty if group is found to be enriched
#' for label only in condition 1} \item{rt1}{Retention time of each
#' isotopologue in condition 1} \item{rt2}{Retention time of each isotopologue
#' in condition 2} \item{relInts1U}{Mean relative intensities of each
#' isotopologue in the samples of condition 1 treated with unlabeled precursor}
#' \item{relInts1L}{Mean relative intensities of each isotopologue in the
#' samples of condition 1 treated with labeled precursor} \item{relInts2U}{Mean
#' relative intensities of each isotopologue in the samples of condition 2
#' treated with unlabeled precursor} \item{relInts2L}{Mean relative intensities
#' of each isotopologue in the samples of condition 2 treated with labeled
#' precursor} \item{p_value}{p-values from Welch's t-tests comparing relative
#' intensities of each isotopologue in labeled samples of condition 1 vs. those
#' of condition 2} \item{sdRelInts1L}{std dev of relative intensity of each
#' isotopologue peak in labeled samples of condition 1} \item{sdRelInts2L}{std
#' dev of relative intensity of each isotopologue peak in labeled samples of
#' condition 2} \item{absInts1L}{mean absolute ion intensity of each
#' isotopologue peak in labeled samples of condition 1} \item{absInts2L}{mean
#' absolute ion intensity of each isotopologue peak in labeled samples of
#' condition 2}
#' @note %% ~~further notes~~
#' @author %% ~~who you are~~ Xiaojing Huang
#' @seealso %% ~~objects to See Also as \code{\link{help}}, ~~~
#' @references %% ~put references to the literature/web site here ~
#' @keywords ~kwd1 ~kwd2
#' @examples
#' 
#' ## Not run: 
#' ## isoDiffReport = getIsoDiffReport(labelsControl, labelsPerturb, condition1 = 'control', condition2 = 'perturb', classesControl = c('C12', 'C12', 'C13', 'C13'), classesPerturb = c('C12', 'C12', 'C13', 'C13'), labeledSamples = 'C13', varEq = FALSE)
#' ##
#' ## This command generates the isoDiffReport for the example given for \code{\link{getIsoLabelReport}}
#' 
#' 
getIsoDiffReport <- function(labelsData1, labelsData2, 
    condition1, condition2, classes1, classes2, labeledSamples, 
    varEq = FALSE, singleSample = FALSE) {
    compounds1 = unlist(labelsData1$compound)
    compounds2 = unlist(labelsData2$compound)
    rt1 = labelsData1$rt
    rt2 = labelsData2$rt
    rel1U = labelsData1$meanRelU
    rel1L = labelsData1$meanRelL
    sdrel1L = labelsData1$sdRelL
    rel2U = labelsData2$meanRelU
    rel2L = labelsData2$meanRelL
    sdrel2L = labelsData2$sdRelL
    abs1L = labelsData1$meanAbsL
    abs2L = labelsData2$meanAbsL
    pattern1 = labelsData1$isotopologue
    pattern2 = labelsData2$isotopologue
    sampData1 = labelsData1$sampleData
    sampData2 = labelsData2$sampleData
    n1 = length(compounds1)
    n2 = length(compounds2)
    n2hasMatched = rep(NA, n2)
    n1_labeled = length(which(classes1 == labeledSamples))
    n2_labeled = length(which(classes2 == labeledSamples))
    base = list()
    class1 = list()
    class2 = list()
    RT1 = list()
    RT2 = list()
    relints1u = list()
    relints1l = list()
    relints2u = list()
    relints2l = list()
    SDrelints1l = list()
    SDrelints2l = list()
    absints1l = list()
    absints2l = list()
    pvalues = list()
    for (i in 1:n1) {
        match1found = 0
        searchCompounds2 = compounds2[is.na(n2hasMatched)]
        len2 = length(searchCompounds2)
        for (j in 1:len2) {
            if (compounds1[i] == searchCompounds2[j]) {
                if (length(pattern1[[i]]) != length(pattern2[is.na(n2hasMatched)][[j]]) || 
                  any(pattern1[[i]] != pattern2[is.na(n2hasMatched)][[j]])) {
                  p1 = pattern1[[i]]
                  p2 = pattern2[is.na(n2hasMatched)][[j]]
                  consensus = unique(c(p1, p2))
                  consensus = consensus[order(consensus)]
                  c1 = list()
                  c2 = list()
                  Ints1 = list()
                  Ints2 = list()
                  rT1 = list()
                  rT2 = list()
                  RI1u = list()
                  RI1l = list()
                  RI2u = list()
                  RI2l = list()
                  sdr1l = list()
                  sdr2l = list()
                  AI1l = list()
                  AI2l = list()
                  for (k in 1:length(consensus)) {
                    P1 = match(consensus[k], p1)
                    P2 = match(consensus[k], p2)
                    if (!is.na(P1)) {
                      Ints1 = c(Ints1, list(sampData1[[i]][P1, 
                        ]))
                      c1 = c(c1, pattern1[[i]][P1])
                      rT1 = c(rT1, rt1[[i]][P1])
                      RI1u = c(RI1u, rel1U[[i]][P1])
                      RI1l = c(RI1l, rel1L[[i]][P1])
                      sdr1l = c(sdr1l, sdrel1L[[i]][P1])
                      AI1l = c(AI1l, abs1L[[i]][P1])
                    } else {
                      Ints1 = c(Ints1, list(rep(0, 
                        dim(sampData1[[i]])[2])))
                      c1 = c(c1, NA)
                      rT1 = c(rT1, NA)
                      RI1u = c(RI1u, NA)
                      RI1l = c(RI1l, NA)
                      sdr1l = c(sdr1l, NA)
                      AI1l = c(AI1l, NA)
                    }
                    if (!is.na(P2)) {
                      Ints2 = c(Ints2, list(sampData2[is.na(n2hasMatched)][[j]][P2, 
                        ]))
                      c2 = c(c2, pattern2[is.na(n2hasMatched)][[j]][P2])
                      rT2 = c(rT2, rt2[is.na(n2hasMatched)][[j]][P2])
                      RI2u = c(RI2u, rel2U[is.na(n2hasMatched)][[j]][P2])
                      RI2l = c(RI2l, rel2L[is.na(n2hasMatched)][[j]][P2])
                      sdr2l = c(sdr2l, sdrel2L[is.na(n2hasMatched)][[j]][P2])
                      AI2l = c(AI2l, abs2L[is.na(n2hasMatched)][[j]][P2])
                    } else {
                      Ints2 = c(Ints2, list(rep(0, 
                        dim(sampData2[is.na(n2hasMatched)][[j]])[2])))
                      c2 = c(c2, NA)
                      rT2 = c(rT2, NA)
                      RI2u = c(RI2u, NA)
                      RI2l = c(RI2l, NA)
                      sdr2l = c(sdr2l, NA)
                      AI2l = c(AI2l, NA)
                    }
                  }
                  class1 = c(class1, list(unlist(c1)))
                  class2 = c(class2, list(unlist(c2)))
                  Ints1 = matrix(unlist(Ints1), nrow = length(consensus), 
                    ncol = dim(sampData1[[i]])[2], 
                    byrow = TRUE)
                  Ints2 = matrix(unlist(Ints2), nrow = length(consensus), 
                    ncol = dim(sampData2[is.na(n2hasMatched)][[j]])[2], 
                    byrow = TRUE)
                  RT1 = c(RT1, list(unlist(rT1)))
                  RT2 = c(RT2, list(unlist(rT2)))
                  relints1u = c(relints1u, list(unlist(RI1u)))
                  relints1l = c(relints1l, list(unlist(RI1l)))
                  relints2u = c(relints2u, list(unlist(RI2u)))
                  relints2l = c(relints2l, list(unlist(RI2l)))
                  SDrelints1l = c(SDrelints1l, list(unlist(sdr1l)))
                  SDrelints2l = c(SDrelints2l, list(unlist(sdr2l)))
                  absints1l = c(absints1l, list(unlist(AI1l)))
                  absints2l = c(absints2l, list(unlist(AI2l)))
                } else {
                  class1 = c(class1, pattern1[i])
                  class2 = c(class2, pattern2[is.na(n2hasMatched)][j])
                  Ints1 = sampData1[[i]]
                  Ints2 = sampData2[is.na(n2hasMatched)][[j]]
                  RT1 = c(RT1, rt1[i])
                  RT2 = c(RT2, rt2[is.na(n2hasMatched)][j])
                  relints1u = c(relints1u, rel1U[i])
                  relints1l = c(relints1l, rel1L[i])
                  relints2u = c(relints2u, rel2U[is.na(n2hasMatched)][j])
                  relints2l = c(relints2l, rel2L[is.na(n2hasMatched)][j])
                  SDrelints1l = c(SDrelints1l, sdrel1L[i])
                  SDrelints2l = c(SDrelints2l, sdrel2L[is.na(n2hasMatched)][j])
                  absints1l = c(absints1l, abs1L[i])
                  absints2l = c(absints2l, abs2L[is.na(n2hasMatched)][j])
                }
                absInts1_labeled = Ints1[, which(classes1 == 
                  labeledSamples), drop = FALSE]
                absInts2_labeled = Ints2[, which(classes2 == 
                  labeledSamples), drop = FALSE]
                relInts1_labeled = absInts1_labeled/matrix(rep(colSums(absInts1_labeled), 
                  dim(Ints1)[1]), nrow = dim(Ints1)[1], 
                  byrow = TRUE)
                relInts2_labeled = absInts2_labeled/matrix(rep(colSums(absInts2_labeled), 
                  dim(Ints2)[1]), nrow = dim(Ints2)[1], 
                  byrow = TRUE)
                if (singleSample == FALSE) {
                  pvalue = list()
                  for (l in 1:dim(relInts1_labeled)[1]) {
                    T = try(t.test(relInts1_labeled[l, 
                      ], relInts2_labeled[l, ], var.equal = varEq), 
                      silent = TRUE)
                    if (class(T) == "try-error" || 
                      is.na(T$p.value)) {
                      pvalue = c(pvalue, 1)
                    } else {
                      pvalue = c(pvalue, T$p.value)
                    }
                  }
                  pvalues = c(pvalues, list(unlist(pvalue)))
                } else {
                  RI1l = unlist(relints1l[[length(relints1l)]])
                  RI2l = unlist(relints2l[[length(relints2l)]])
                  RI1l[which(is.na(RI1l))] = 0
                  RI2l[which(is.na(RI2l))] = 0
                  pvalue = sum(abs(RI1l - RI2l))
                  pvalues = c(pvalues, pvalue)
                }
                base = c(base, compounds1[i])
                n2hasMatched[is.na(n2hasMatched)][j] = 1
                match1found = 1
                break
            }
        }
        if (!match1found) {
            base = c(base, compounds1[i])
            class1 = c(class1, pattern1[i])
            class2 = c(class2, -1)
            RT1 = c(RT1, rt1[i])
            RT2 = c(RT2, -1)
            pvalues = c(pvalues, 0)
            relints1u = c(relints1u, rel1U[i])
            relints1l = c(relints1l, rel1L[i])
            relints2u = c(relints2u, -1)
            relints2l = c(relints2l, -1)
            SDrelints1l = c(SDrelints1l, sdrel1L[i])
            SDrelints2l = c(SDrelints2l, -1)
            absints1l = c(absints1l, abs1L[i])
            absints2l = c(absints2l, -1)
        }
    }
    unmatchedCompounds2 = compounds2[is.na(n2hasMatched)]
    len2 = length(unmatchedCompounds2)
    for (i in 1:len2) {
        base = c(base, unmatchedCompounds2[i])
        class1 = c(class1, -1)
        class2 = c(class2, pattern2[is.na(n2hasMatched)][i])
        RT1 = c(RT1, -1)
        RT2 = c(RT2, rt2[is.na(n2hasMatched)][i])
        pvalues = c(pvalues, 0)
        relints1u = c(relints1u, -1)
        relints1l = c(relints1l, -1)
        relints2u = c(relints2u, rel2U[is.na(n2hasMatched)][i])
        relints2l = c(relints2l, rel2L[is.na(n2hasMatched)][i])
        SDrelints1l = c(SDrelints1l, -1)
        SDrelints2l = c(SDrelints2l, sdrel2L[is.na(n2hasMatched)][i])
        absints1l = c(absints1l, -1)
        absints2l = c(absints2l, abs2L[is.na(n2hasMatched)][i])
    }
    isoDiffReport = list(class1[order(unlist(base))], 
        class2[order(unlist(base))], RT1[order(unlist(base))], 
        RT2[order(unlist(base))], relints1u[order(unlist(base))], 
        relints1l[order(unlist(base))], relints2u[order(unlist(base))], 
        relints2l[order(unlist(base))], pvalues[order(unlist(base))], 
        SDrelints1l[order(unlist(base))], SDrelints2l[order(unlist(base))], 
        absints1l[order(unlist(base))], absints2l[order(unlist(base))])
    names(isoDiffReport) = c(condition1, condition2, 
        "rt1", "rt2", "relInts1U", "relInts1L", "relInts2U", 
        "relInts2L", "p_value", "sdrelInts1L", "sdrelInts2L", 
        "absInts1L", "absInts2L")
    return(isoDiffReport)
}


#' %% ~~function to do ... ~~ Print X13CMS reports.
#' 
#' %% ~~ A concise (1-5 lines) description of what the function does. ~~ Prints
#' either an isoLabelReport or isoDiffReport to tab-delimited text file. The
#' column headings are the names of the list components in the report.
#' 
#' %% ~~ If necessary, more details than the description above ~~
#' 
#' @param listReport %% ~~Describe \code{listReport} here~~ either an
#' isoLabelReport or isoDiffReport
#' @param outputfile %% ~~Describe \code{outputfile} here~~ character variable
#' of output file's name
#' @return %% ~Describe the value returned %% If it is a LIST, use %%
#' \item{comp1 }{Description of 'comp1'} %% \item{comp2 }{Description of
#' 'comp2'} %% ... None. Side effect only.
#' @note %% ~~further notes~~
#' @author %% ~~who you are~~ Xiaojing Huang
#' @seealso %% ~~objects to See Also as \code{\link{help}}, ~~~
#' @references %% ~put references to the literature/web site here ~
#' @keywords ~kwd1 ~kwd2
#' @examples
#' 
#' 
#' 
printIsoListOutputs <- function(listReport, outputfile) {
    colNames = names(listReport)
    nblocks = length(colNames)
    colNames = c("", colNames)
    write(colNames, file = outputfile, ncolumns = nblocks + 
        1, append = FALSE, sep = "\t")
    nrows = length(listReport[[1]])
    for (i in 1:nrows) {
        maxRows = 1
        ncols = 0
        for (j in 1:nblocks) {
            cell = listReport[[j]][[i]]
            if (length(cell) > maxRows && all(is.numeric(cell))) {
                maxRows = length(cell)
                ncols = ncols + 1
            } else if (all(is.matrix(cell))) {
                ncols = ncols + dim(cell)[2]
            } else {
                ncols = ncols + 1
            }
        }
        rowMatrix = matrix(rep(NA, (maxRows + 1) * 
            ncols), nrow = maxRows + 1)
        for (j in 1:maxRows) {
            ncol = 1
            for (k in 1:nblocks) {
                cell = listReport[[k]][[i]]
                if (length(cell) == 1 && j == 1 && 
                  all(cell != -1)) {
                  rowMatrix[j, ncol] = cell
                  ncol = ncol + 1
                } else if (all(is.numeric(cell)) || 
                  all(is.integer(cell)) || all(is.character(cell))) {
                  if (!all(is.na(cell[j])) && all(cell[j] != 
                    -1)) {
                    rowMatrix[j, ncol] = cell[j]
                  }
                  ncol = ncol + 1
                } else if (all(is.matrix(cell))) {
                  rowMatrix[j, ncol:(ncol + dim(cell)[2] - 
                    1)] = cell[j, ]
                  ncol = ncol + dim(cell)[2]
                }
            }
        }
        nr = dim(rowMatrix)[1]
        rowMatrix = cbind(c(i, rep(NA, nr - 1)), rowMatrix)
        write.table(rowMatrix, outputfile, na = "", 
            row.names = FALSE, col.names = FALSE, 
            sep = "\t", append = TRUE)
    }
}


#' %% ~~function to do ... ~~ Plot isotopologue groups from labeling report
#' 
#' %% ~~ A concise (1-5 lines) description of what the function does. ~~
#' Generates bar plots of ion intensity distribution among the peaks of an
#' isotopologue group
#' 
#' %% ~~ If necessary, more details than the description above ~~
#' 
#' @param isoLabelReport %% ~~Describe \code{isoLabelReport} here~~ Output of
#' getIsoLabelReport().
#' @param intOption %% ~~Describe \code{intOption} here~~ Choice of plotting
#' relative ('rel') or absolute ('abs') ion intensities. Defaults to 'rel'.
#' @param classes Character vector designating whether each sample is unlabeled
#' (e.g. 'C12') or labeled (e.g. 'C13')
#' @param labeledSamples %% ~~Describe \code{labeledSamples} here~~ character
#' variable designating labeled samples (e.g. 'C13')
#' @param outputfile %% ~~Describe \code{outputfile} here~~ Name of pdf file to
#' which plots are drawn.
#' @return %% ~Describe the value returned %% If it is a LIST, use %%
#' \item{comp1 }{Description of 'comp1'} %% \item{comp2 }{Description of
#' 'comp2'} %% ... Pdf file containing 12 isotopologue group plots per page.
#' @note %% ~~further notes~~
#' @author %% ~~who you are~~ Xiaojing Huang
#' @seealso %% ~~objects to See Also as \code{\link{help}}, ~~~
#' @references %% ~put references to the literature/web site here ~
#' @keywords ~kwd1 ~kwd2
#' @examples
#' 
#' ## Not run:
#' ## plotLabelReport(labelsControl, 'rel', 'labelsControlPlots.pdf')
#' 
plotLabelReport <- function(isoLabelReport, intOption = "rel", 
    classes, labeledSamples, outputfile) {
    pdf(outputfile, width = 8.5, height = 11)
    numGroups = length(isoLabelReport[[1]])
    par(mfrow = c(4, 3))
    for (i in 1:numGroups) {
        if (i%%12 == 1) {
            par(mfrow = c(4, 3))
        }
        if (intOption == "rel") {
            toPlot = cbind(isoLabelReport$meanRelU[[i]], 
                isoLabelReport$meanRelL[[i]])
            sds = cbind(isoLabelReport$sdRelU[[i]], 
                isoLabelReport$sdRelL[[i]])
            YLIM = c(0, 1)
        } else if (intOption == "abs") {
            toPlot = cbind(isoLabelReport$meanAbsU[[i]], 
                isoLabelReport$meanAbsL[[i]])
            sds = cbind(isoLabelReport$cvTotalU[[i]] * 
                isoLabelReport$meanAbsU[[i]], isoLabelReport$cvTotalL[[i]] * 
                isoLabelReport$meanAbsL[[i]])
            YLIM = NULL
            totalPoolsU = colSums(isoLabelReport$sampleData[[i]][, 
                which(classes != labeledSamples)])
            totalPoolsL = colSums(isoLabelReport$sampleData[[i]][, 
                which(classes == labeledSamples)])
            T = t.test(totalPoolsU, totalPoolsL)
            pval = round(T$p.value, 3)
        } else {
            print("Error. Specify intOption as rel or abs.")
            break
        }
        numPeaks = dim(toPlot)[1]
        rownames(toPlot) = as.character(round(isoLabelReport$isotopologue[[i]], 
            4))
        colnames(toPlot) = c("U", "L")
        if (intOption == "rel") {
            mp = barplot(t(toPlot), beside = TRUE, 
                legend.text = c("U", "L"), main = paste("m/z = ", 
                  rownames(toPlot)[1], "\n", "RT = ", 
                  isoLabelReport$rt[[i]][1]), axisnames = FALSE, 
                ylim = YLIM)
            text(seq(1.5, (numPeaks - 1) * 3 + 1.5, 
                by = 3), par("usr")[3] - 0.2, xpd = TRUE, 
                labels = rownames(toPlot), srt = 45)
            segments(mp, t(toPlot) - t(sds), mp, t(toPlot) + 
                t(sds), lwd = 2)
        } else {
            barplot(toPlot, legend.text = rownames(toPlot), 
                main = paste("m/z = ", rownames(toPlot)[1], 
                  "\n", "RT = ", isoLabelReport$rt[[i]][1], 
                  "\n", "p = ", pval))
        }
    }
    dev.off()
}


#' %% ~~function to do ... ~~ Plot output of getIsoDiffReport()
#' 
#' %% ~~ A concise (1-5 lines) description of what the function does. ~~
#' Generates bar plots comparing isotopologue distributions in labeled samples
#' of two different sample conditions. Use relative ion intensities to show
#' distributions.
#' 
#' %% ~~ If necessary, more details than the description above ~~
#' 
#' @param isoDiffReport %% ~~Describe \code{isoDiffReport} here~~ Output of
#' getIsoDiffReport().
#' @param xcmsSet xcmsSet object containing grouped and retention-time-aligned
#' peaks (i.e. after calling group() and retcor() in XCMS)
#' @param intChoice 'intChoice' from getIsoLabelReport()
#' @param sampleNames1 'sampleNames' from getIsoLabelReport() called for
#' condition 1
#' @param sampleNames2 'sampleNames' from getIsoLabelReport() called for
#' condition 2
#' @param labelReport1 output of getIsoLabelReport() for condition 1
#' @param labelReport2 output of getIsoLabelReport() for condition 2
#' @param classes1 character vector designating whether each sample in
#' labelsReport1 is unlabeled or labeled
#' @param classes2 character vector designating whether each sample in
#' labelsReport2 is unlabeled or labeled
#' @param labeledSamples character variable designating labeled samples (e.g.
#' 'C13')
#' @param isotopeMassDifference mass difference between unlabeled and labeled
#' atom, e.g. 1.00335 for C13
#' @param outputfile %% ~~Describe \code{outputfile} here~~ Name of pdf file to
#' which plots are drawn.
#' @return %% ~Describe the value returned %% If it is a LIST, use %%
#' \item{comp1 }{Description of 'comp1'} %% \item{comp2 }{Description of
#' 'comp2'} %% ... Pdf file containing 12 comparison plots per page.
#' @note %% ~~further notes~~
#' @author %% ~~who you are~~ Xiaojing Huang
#' @seealso %% ~~objects to See Also as \code{\link{help}}, ~~~
#' @references %% ~put references to the literature/web site here ~
#' @keywords ~kwd1 ~kwd2
#' @examples
#' 
#' ## Not run:
#' ## plotIsoDiffReport(isoDiffReport, 'isoDiffPlots.pdf') 
#' 
plotIsoDiffReport <- function(isoDiffReport, xcmsSet, 
    intChoice, sampleNames1, sampleNames2, labelReport1, 
    labelReport2, classes1, classes2, labeledSamples, 
    isotopeMassDifference, outputfile) {
    pdf(outputfile, width = 8.5, height = 11)
    numGroups = length(isoDiffReport[[1]])
    par(mfrow = c(4, 3))
    ints = groupval(xcmsSet, "medret", intChoice)
    ints1 = ints[, which(colnames(ints) %in% sampleNames1[which(classes1 == 
        labeledSamples)])]
    ints2 = ints[, which(colnames(ints) %in% sampleNames2[which(classes2 == 
        labeledSamples)])]
    ints1[is.na(ints1)] = 0
    ints2[is.na(ints2)] = 0
    for (i in 1:numGroups) {
        if (i%%12 == 1) {
            par(mfrow = c(4, 3))
        }
        if (isoDiffReport$relInts1L[[i]][1] != -1 & 
            isoDiffReport$relInts2L[[i]][1] != -1) {
            toPlot = cbind(isoDiffReport$relInts1L[[i]], 
                isoDiffReport$relInts2L[[i]])
            sds = cbind(isoDiffReport$sdrelInts1L[[i]], 
                isoDiffReport$sdrelInts2L[[i]])
            numPeaks = dim(toPlot)[1]
            a = isoDiffReport[[1]][[i]]
            b = isoDiffReport[[2]][[i]]
            if (NA %in% a) {
                mz = round(b[1], 4)
                isos = round((b - b[1])/isotopeMassDifference)
                rownames(toPlot) = as.character(isos)
            } else {
                mz = round(a[1], 4)
                isos = round((a - a[1])/isotopeMassDifference)
                rownames(toPlot) = as.character(isos)
            }
            mp = barplot(t(toPlot), beside = TRUE, 
                legend.text = names(isoDiffReport)[1:2], 
                main = paste("m/z = ", mz, "\n", "RT = ", 
                  isoDiffReport$rt1[[i]][1]))
            sds[which(is.na(sds))] = 0
            segments(mp, t(toPlot) - t(sds), mp, t(toPlot) + 
                t(sds), lwd = 2)
        } else if (isoDiffReport$relInts1L[[i]][1] == 
            -1) {
            a = length(isoDiffReport$relInts2L[[i]])
            b = isoDiffReport[[2]][[i]]
            ind = which(unlist(labelReport2$compound) == 
                b[1])
            abs = ints1[labelReport2$groupID[[ind]], 
                ]
            tots = colSums(abs)
            rel = abs/matrix(rep(tots, a), byrow = TRUE, 
                nrow = a)
            meanRel = rowMeans(rel)
            sd = apply(rel, 1, sd)
            mz = round(isoDiffReport[[2]][[i]][1], 
                4)
            toPlot = cbind(meanRel, isoDiffReport$relInts2L[[i]])
            sds = cbind(sd, isoDiffReport$sdrelInts2L[[i]])
            numPeaks = dim(toPlot)[1]
            isos = round((isoDiffReport[[2]][[i]] - 
                isoDiffReport[[2]][[i]][1])/isotopeMassDifference)
            rownames(toPlot) = as.character(isos)
            mp = barplot(t(toPlot), beside = TRUE, 
                legend.text = names(isoDiffReport)[1:2], 
                main = paste("m/z = ", mz, "\n", "RT = ", 
                  isoDiffReport$rt2[[i]][1]))
            sds[which(is.na(sds))] = 0
            segments(mp, t(toPlot) - t(sds), mp, t(toPlot) + 
                t(sds), lwd = 2)
        } else {
            a = length(isoDiffReport$relInts1L[[i]])
            b = isoDiffReport[[1]][[i]]
            ind = which(unlist(labelReport1$compound) == 
                b[1])
            abs = ints2[labelReport1$groupID[[ind]], 
                ]
            tots = colSums(abs)
            rel = abs/matrix(rep(tots, a), byrow = TRUE, 
                nrow = a)
            meanRel = rowMeans(rel)
            sd = apply(rel, 1, sd)
            mz = round(isoDiffReport[[1]][[i]][1], 
                4)
            toPlot = cbind(isoDiffReport$relInts1L[[i]], 
                meanRel)
            sds = cbind(isoDiffReport$sdrelInts1L[[i]], 
                sd)
            numPeaks = dim(toPlot)[1]
            isos = round((isoDiffReport[[1]][[i]] - 
                isoDiffReport[[1]][[i]][1])/isotopeMassDifference)
            rownames(toPlot) = as.character(isos)
            mp = barplot(t(toPlot), beside = TRUE, 
                legend.text = names(isoDiffReport)[1:2], 
                main = paste("m/z = ", mz, "\n", "RT = ", 
                  isoDiffReport$rt1[[i]][1]))
            sds[which(is.na(sds))] = 0
            segments(mp, t(toPlot) - t(sds), mp, t(toPlot) + 
                t(sds), lwd = 2)
        }
    }
    dev.off()
}


#' %% ~~function to do ... ~~ Plot output of getIsoDiffReport()
#' 
#' %% ~~ A concise (1-5 lines) description of what the function does. ~~
#' Generates bar plots comparing isotopologue distributions in labeled samples
#' of two different sample classes. Use absolute ion intensities to show
#' distributions and compares total pool sizes of isotopologue group between
#' sample classes.
#' 
#' %% ~~ If necessary, more details than the description above ~~
#' 
#' @param isoDiffReport %% ~~Describe \code{isoDiffReport} here~~ Output of
#' getIsoDiffReport().
#' @param xcmsSet xcmsSet object containing grouped and retention-time-aligned
#' peaks (i.e. after calling group() and retcor() in XCMS)
#' @param intChoice 'intChoice' from getIsoLabelReport()
#' @param sampleNames1 'sampleNames' from getIsoLabelReport() called for
#' condition 1
#' @param sampleNames2 'sampleNames' from getIsoLabelReport() called for
#' condition 2
#' @param labelReport1 output of getIsoLabelReport() for condition 1
#' @param labelReport2 output of getIsoLabelReport() for condition 2
#' @param classes1 character vector designating whether each sample in
#' labelsReport1 is unlabeled or labeled
#' @param classes2 character vector designating whether each sample in
#' labelsReport2 is unlabeled or labeled
#' @param labeledSamples character variable designating labeled samples (e.g.
#' 'C13')
#' @param outputfile %% ~~Describe \code{outputfile} here~~ Name of pdf file to
#' which plots are drawn.
#' @return %% ~Describe the value returned %% If it is a LIST, use %%
#' \item{comp1 }{Description of 'comp1'} %% \item{comp2 }{Description of
#' 'comp2'} %% ... Pdf file containing 12 comparison plots per page.
#' @note %% ~~further notes~~
#' @author %% ~~who you are~~ Xiaojing Huang
#' @seealso %% ~~objects to See Also as \code{\link{help}}, ~~~
#' @references %% ~put references to the literature/web site here ~
#' @keywords ~kwd1 ~kwd2
#' @examples
#' 
#' ## Not run:
#' ## plotTotalPools(isoDiffReport, 'totalPools.pdf') 
#' 
plotTotalIsoPools <- function(isoDiffReport, xcmsSet, 
    intChoice, sampleNames1, sampleNames2, labelReport1, 
    labelReport2, classes1, classes2, labeledSamples, 
    outputfile) {
    pdf(outputfile, width = 8.5, height = 11)
    numGroups = length(isoDiffReport[[1]])
    par(mfrow = c(4, 3))
    ints = groupval(xcmsSet, "medret", intChoice)
    ints1 = ints[, which(colnames(ints) %in% sampleNames1[which(classes1 == 
        labeledSamples)])]
    ints2 = ints[, which(colnames(ints) %in% sampleNames2[which(classes2 == 
        labeledSamples)])]
    ints1[is.na(ints1)] = 0
    ints2[is.na(ints2)] = 0
    for (i in 1:numGroups) {
        if (i%%12 == 1) {
            par(mfrow = c(4, 3))
        }
        if (isoDiffReport$absInts1L[[i]][1] != -1 & 
            isoDiffReport$absInts2L[[i]][1] != -1) {
            toPlot = cbind(isoDiffReport$absInts1L[[i]], 
                isoDiffReport$absInts2L[[i]])
            toPlot[is.na(toPlot)] = 0
            a = isoDiffReport[[1]][[i]]
            b = isoDiffReport[[2]][[i]]
            ind1 = which(unlist(labelReport1$compound) == 
                a[1])
            ind2 = which(unlist(labelReport2$compound) == 
                b[1])
            totalPools1 = colSums(labelReport1$sampleData[[ind1]][, 
                which(classes1 == labeledSamples)])
            totalPools2 = colSums(labelReport2$sampleData[[ind2]][, 
                which(classes2 == labeledSamples)])
            T = t.test(totalPools1, totalPools2)
            pval = round(T$p.value, 3)
            if (NA %in% a) {
                rownames(toPlot) = as.character(round(b - 
                  b[1]))
            } else {
                rownames(toPlot) = as.character(round(a - 
                  a[1]))
            }
            colnames(toPlot) = c(unlist(attributes(isoDiffReport))[1], 
                unlist(attributes(isoDiffReport))[2])
            barplot(toPlot, legend.text = rownames(toPlot), 
                main = paste("m/z =", round(b[1], 
                  4), ", ", "RT =", isoDiffReport$rt1[[i]][1], 
                  "\n", "p =", pval))
        } else if (isoDiffReport$absInts1L[[i]][1] == 
            -1) {
            b = isoDiffReport[[2]][[i]]
            ind = which(unlist(labelReport2$compound) == 
                b[1])
            a = rowMeans(ints1[labelReport2$groupID[[ind]], 
                ])
            totalPools1 = colSums(ints1[labelReport2$groupID[[ind]], 
                ])
            totalPools2 = colSums(labelReport2$sampleData[[ind]][, 
                which(classes2 == labeledSamples)])
            T = t.test(totalPools1, totalPools2)
            pval = round(T$p.value, 3)
            toPlot = cbind(a, isoDiffReport$absInts2L[[i]])
            rownames(toPlot) = as.character(round(isoDiffReport[[2]][[i]] - 
                isoDiffReport[[2]][[i]][1]))
            colnames(toPlot) = c(unlist(attributes(isoDiffReport))[1], 
                unlist(attributes(isoDiffReport))[2])
            barplot(toPlot, legend.text = rownames(toPlot), 
                main = paste("m/z =", round(isoDiffReport[[2]][[i]][1], 
                  4), ", ", "RT =", isoDiffReport$rt2[[i]][1], 
                  "\n", "p =", pval, "\n", "(", colnames(toPlot)[1], 
                  "not significantly labeled )"))
        } else {
            b = isoDiffReport[[1]][[i]]
            ind = which(unlist(labelReport1$compound) == 
                b[1])
            a = rowMeans(ints2[labelReport1$groupID[[ind]], 
                ])
            totalPools2 = colSums(ints2[labelReport1$groupID[[ind]], 
                ])
            totalPools1 = colSums(labelReport1$sampleData[[ind]][, 
                which(classes1 == labeledSamples)])
            T = t.test(totalPools1, totalPools2)
            pval = round(T$p.value, 3)
            toPlot = cbind(isoDiffReport$absInts1L[[i]], 
                a)
            rownames(toPlot) = as.character(round(isoDiffReport[[1]][[i]] - 
                isoDiffReport[[1]][[i]][1]))
            colnames(toPlot) = c(unlist(attributes(isoDiffReport))[1], 
                unlist(attributes(isoDiffReport))[2])
            barplot(toPlot, legend.text = rownames(toPlot), 
                main = paste("m/z =", round(isoDiffReport[[1]][[i]][1], 
                  4), ", ", "RT =", isoDiffReport$rt1[[i]][1], 
                  "\n", "p =", pval, "\n", "(", colnames(toPlot)[2], 
                  "not significantly labeled )"))
        }
    }
    dev.off()
}


#' %% ~~function to do ... ~~ Get XCMS-style diffReport from X13CMS data.
#' 
#' %% ~~ A concise (1-5 lines) description of what the function does. ~~
#' Generates conventional XCMS-like diff report comparing peaks between two
#' sample classes, e.g. 'control', and 'treatment' from an xcmsSet object
#' configured for X13CMS (i.e. one in which the sample classes have to be
#' demarcated as unlabeled and labeled)
#' 
#' %% ~~ If necessary, more details than the description above ~~
#' 
#' @param xcmsSet %% ~~Describe \code{xcmsSet} here~~ XCMS object
#' @param class1sampNames %% ~~Describe \code{class1sampNames} here~~ The
#' subset of rownames(xcmsSet@phenoData) corresponding to names of the
#' unlabeled samples of class 1 (e.g. 'control').
#' @param class2sampNames %% ~~Describe \code{class2sampNames} here~~ The
#' subset of rownames(xcmsSet@phenoData) corresponding to names of the
#' unlabeled samples of class 2 (e.g. 'perturb').
#' @param varEq %% ~~Describe \code{varEq} here~~ Boolean indicating whether to
#' assume that absolute ion intensities of each peak in both sample classes are
#' drawn from distributions with equal variance. Defaults to FALSE.
#' @param intChoice %% ~~Describe \code{intChoice} here~~ one of 'maxo',
#' 'into', or 'intb'--the choice of which peak intensity measurement to use
#' from the XCMS object.
#' @return %% ~Describe the value returned %% If it is a LIST, use %%
#' \item{comp1 }{Description of 'comp1'} %% \item{comp2 }{Description of
#' 'comp2'} A dataframe whose rows correspond to peak groups identified in the
#' xcmsSet object and whose columns are: - pvalue: p-value of Welch's t-test
#' comparing absolute intensities of the peak in the two conditions -
#' foldChange: change in mean peak intensity in condition 2 compared to
#' conditions 1 - mzmed: median m/z of peak - rtmed: median retention time of
#' peak - means1: mean peak intensity in condition 1 - means2: mean peak
#' intensity in condition 2 - peakIntensities1: all peak intensities in the
#' group in condition 1 - peakIntensities2: all peak intensities in the group
#' in condition 2 %% ...
#' @note %% ~~further notes~~
#' @author %% ~~who you are~~ Xiaojing Huang
#' @seealso %% ~~objects to See Also as \code{\link{help}}, ~~~
#' @references %% ~put references to the literature/web site here ~
#' @keywords ~kwd1 ~kwd2
#' @examples
#' 
#' ## Not run:
#' 
#' ## diffReport = miniDiffReport(xcmsSet, class1sampNames, class2sampNames, varEq = FALSE, intChoice = 'intb') 
#' 
#' ## From the example given in getIsoLabelReport(), the variables class1sampNames and class2sampNames would be c('control_unlabeled1', 'control_unlabeled2'), and c('perturb_unlabeled1', 'perturb_unlabeled2'), respectively.
#' 
miniDiffReport <- function(xcmsSet, class1sampNames, 
    class2sampNames, varEq = FALSE, intChoice) {
    peakIntensities = groupval(xcmsSet, method = "medret", 
        value = intChoice)
    if (intChoice == "intb") {
        peakIntensities[is.na(peakIntensities)] = 0
    }
    peakIntensities1 = peakIntensities[, match(class1sampNames, 
        colnames(peakIntensities))]
    peakIntensities2 = peakIntensities[, match(class2sampNames, 
        colnames(peakIntensities))]
    means1 = rowMeans(peakIntensities1)
    means2 = rowMeans(peakIntensities2)
    groups = data.frame(xcmsSet@groups)
    npeaks = dim(peakIntensities)[1]
    pvalue = rep(0, npeaks)
    foldChange = rep(0, npeaks)
    for (i in 1:npeaks) {
        T = t.test(peakIntensities1[i, ], peakIntensities2[i, 
            ], var.equal = varEq)
        pvalue[i] = T$p.value
        if (means2[i] > means1[i]) {
            foldChange[i] = means2[i]/means1[i]
        } else {
            foldChange[i] = -means1[i]/means2[i]
        }
    }
    foldChange = foldChange[order(pvalue)]
    peakIntensities1 = peakIntensities1[order(pvalue), 
        ]
    peakIntensities2 = peakIntensities2[order(pvalue), 
        ]
    means1 = means1[order(pvalue)]
    means2 = means2[order(pvalue)]
    groups = groups[order(pvalue), c(1, 4)]
    pvalue = pvalue[order(pvalue)]
    out = cbind(pvalue, foldChange, groups, means1, 
        means2, peakIntensities1, peakIntensities2)
    return(out)
}


#' %% ~~function to do ... ~~ Filter isoDiffReport
#' 
#' %% ~~ A concise (1-5 lines) description of what the function does. ~~ This
#' routine filters an isoDiffReport to include only isotopologue groups that
#' are differentially labeled between sample conditions and are not falsely
#' unique due to suboptimal chromatography conditions. An additional, optional
#' filter exists for the special case of 13C labeling experiments; it removes
#' groups that consist only of M and M+1 peaks.
#' 
#' %% ~~ If necessary, more details than the description above ~~
#' 
#' @param isoDiffReport %% ~~Describe \code{isoDiffReport} here~~ report to be
#' filtered
#' @param alpha %% ~~Describe \code{alpha} here~~ p-value cutoff to call
#' labeling pattern of isotoplogue group as different between sample
#' conditions.
#' @param whichPeak %% ~~Describe \code{whichPeak} here~~ option for which
#' isotopologue peaks within the group to consider when calling significance
#' for differential labeling. Default = 1 indicates consideration of the
#' relative intensities of the base peak only. '2' indicates consideration of
#' all peaks and a call of significance if at least one of their relative
#' intensities is significantly different between the two sample conditions.
#' @return %% ~Describe the value returned %% If it is a LIST, use %%
#' \item{comp1 }{Description of 'comp1'} %% \item{comp2 }{Description of
#' 'comp2'} %% ... isoDiffReport
#' @note %% ~~further notes~~
#' @author %% ~~who you are~~ Xiaojing Huang
#' @seealso %% ~~objects to See Also as \code{\link{help}}, ~~~
#' @references %% ~put references to the literature/web site here ~
#' @keywords ~kwd1 ~kwd2
#' @examples
#' 
#' ## Not run:
#' ## filterIsoDiffReport(isoDiffReport, alpha = 0.05, whichPeak = 1, is13C = FALSE)   
#'  
#' 
filterIsoDiffReport <- function(isoDiffReport, alpha, 
    whichPeak = 2) {
    n = length(isoDiffReport[[1]])
    outtakes = list()
    for (i in 1:n) {
        if (any(unlist(isoDiffReport$p_value[i]) == 
            1)) {
            outtakes = c(outtakes, i)
            next
        }
        if (whichPeak == 1 & all(isoDiffReport$p_value[i][[1]] > 
            alpha)) {
            outtakes = c(outtakes, i)
            next
        }
        if (whichPeak == 2 & !any(unlist(isoDiffReport$p_value[i]) < 
            alpha)) {
            outtakes = c(outtakes, i)
            next
        }
    }
    outtakes = unlist(outtakes)
    if (length(outtakes) == 0) {
        out = isoDiffReport
    } else {
        out = lapply(isoDiffReport, "[", -outtakes)
    }
    return(out)
}
