require(xcms)
require(stats)

# extracts table of peak intensities for samples named in "sampleNames"
getPeakIntensities <- function(xcmsSet, sampleNames, intChoice) {
  # intChoice refers to choice of "maxo", "into", or "intb" for peak intensities
  peakIntensities = groupval(xcmsSet, method = "medret", value = intChoice); 
  peakIntensities = peakIntensities[, match(sampleNames, colnames(peakIntensities))];
  if (intChoice == "intb") {
    peakIntensities[is.na(peakIntensities)] = 0;
  }
  return(peakIntensities);
}

# outputs report of features enriched for labeled isotope in samples named in "sampleNames"
getIsoLabelReport <- function(xcmsSet, sampleNames, unlabeledSamples, labeledSamples, isotopeMassDiff, RTwindow, ppm, 
                              massOfLabeledAtom, noiseCutoff, intChoice = "intb", varEq = FALSE, alpha, singleSample = FALSE, compareOnlyDistros = FALSE,
                              monotonicityTol = FALSE, enrichTol = 0.1) {
  # Variables:
  # xcmsSet: xcmsSet object containing grouped and retention-time-aligned peaks (i.e. after calling group() and retcor() in XCMS); 
  #          sample classes should be designated as "unlabeledSamples" and "labeledSamples";
  #          may contain samples from one or more biological conditions
  # sampleNames: vector of names of the unlabeled and labeled samples in the condition for which the report is being generated;
  #              taken from rownames of xcmsSet@phenoData
  # unlabeledSamples: character variable name of unlabeled samples as designated in xcmsSet@phenoData (e.g. "C12")
  # labeledSamples: character variable name of labeled samples (e.g. "C13")
  # isotopeMassDiff: difference in mass between labeled and unlabeled atom (e.g. 1.00335 for C13) 
  # RTwindow: retention time window in which all peaks are considered to be co-eluting
  # ppm: ppm allowance for deviation of peaks within an isotopologue group from expected m/z; instrument dependent
  # massOfLabeledAtom: e.g. 12.0000 for C12
  # noiseCutoff: ion intensity cutoff below which a peak is considered noise; instrument dependent
  # intChoice: choice of "maxo", "into", or "intb" for peak intensities
  # varEq: Boolean indicating whether to assume equal variances for the peak intensities used in the t-test comparing unlabeled
  #        vs labeled samples
  # alpha: p-value cutoff for significance of label enrichment
  # singleSample: Boolean indicating whether only a single sample exists for the unlabeled and labeled conditions; 
  #               if TRUE, all potential isotopologue groups are returned and no statistics are performed
  # compareOnlyDistros: Boolean indicating whether to assess for true labeling patterns (FALSE) or 
  #                     to compare label distributions if both samples are labeled (TRUE)
  # monotonicityTol: tolerance parameter used to enforce expected ion intensity pattern (i.e. monotonic decrease from M0 to Mn) 
  #                  in unlabeled samples; a low value closer to 0 enforces stricter monotonicity; default is to not enforce monotonicity
  #                  due to potential carryover between samples
  # enrichTol: tolerance parameter for enforcing enrichment of higher isotopologues in labeled samples; a value of 0 enforces strict requirement
  #            for enrichment of higher isotopologues to be higher in labeled samples
  
  ### UNPACK THE XCMS OBJECT ###
  groups = data.frame(xcmsSet@groups);
  peakIntensities = getPeakIntensities(xcmsSet, sampleNames, intChoice);
  peakIntensities = peakIntensities[order(groups$mzmed),];
  groups = groups[order(groups$mzmed),]; # order peaks by m/z
  groupRTs = groups$rtmed;
  groupMzs = groups$mzmed;
  groupIDs = as.numeric(rownames(groups));
  nGroups = length(groupMzs);
  classes = as.character(xcmsSet@phenoData[match(sampleNames,rownames(xcmsSet@phenoData)),]);
  numSamples = length(classes);
  
  intensities1= peakIntensities[,which(classes == unlabeledSamples), drop = FALSE];
  intensities2 = peakIntensities[,which(classes == labeledSamples), drop = FALSE];
  
  iMD = isotopeMassDiff; 
  
  ### ISOTOPOLOGUE SEARCH: ###
  # temporary storage lists for potential isotopologue pairs
  base = list(); # index of base peak
  labeled = list(); # index of labeled peak
  basePeak = list(); # m/z of base peak
  labeledPeak = list(); # m/z of labeled peak
  
  groupIndicesByRT = order(groupRTs);
  orderedGroupRTs = groupRTs[groupIndicesByRT];
  # the search:
  for (i in 1:nGroups) {
    binI = groupIndicesByRT[orderedGroupRTs - orderedGroupRTs[i] >= 0 & orderedGroupRTs - orderedGroupRTs[i] <= RTwindow];
    bin = groups[binI,];
    binSize = length(binI);
    I = groupIndicesByRT[i];
    if (binSize > 0) {
      # do pairwise comparisons between every m/z group represented in that bin
      # a = unlabeled peak
      # b = labeled peak
      for (j in 1:binSize) {
        if (groups$mzmed[I] < bin$mzmed[j]) {
          a = I;
          b = binI[j];
        }
        else {
          a = binI[j];
          b = I;
        }     
   
        delta = (groupMzs[b] - groupMzs[a])/iMD;
        DELTA = round(delta); # candidate number of labeled atoms
        if (DELTA == 0) {
          next;
        }
        # if delta is a multiple of mass defect to within ppm error:
        if ( delta <= DELTA*(1+ppm/1000000) + (groupMzs[a]*ppm/1000000)/(iMD*(1-ppm/1000000)) 
             && delta >= DELTA*(1-ppm/1000000) - (groupMzs[a]*ppm/1000000)/(iMD*(1+ppm/1000000)) ) {
          # ...check if mass defect is too large
          if ( DELTA  * massOfLabeledAtom >= groupMzs[a]) {
            next;
          }
          # ...check that intensities of labeled peak are less than those of base peak in the unlabeled samples (class == 1)
          if ( mean(intensities1[b,]) > mean(intensities1[a,]) && !compareOnlyDistros ) {
            next;
          } 
          # ...check that intensities are not 0 in all samples
          if ( all(intensities1[a,] == 0) && all(intensities2[a,] == 0)) {
            next;
          }
          if ( all(intensities1[b,] == 0) && all(intensities2[b,] == 0)) {
            next;
          }
          # record pair of peaks if all filters are passed
          base = c(base, a);
          labeled = c(labeled, b);
          basePeak = c(basePeak, groupMzs[a]);
          labeledPeak = c(labeledPeak, groupMzs[b]);
        }
      }
    }    
  }
  
  labelsMatrix = as.matrix(cbind(unlist(base), 
                                 unlist(labeled), 
                                 unlist(basePeak), 
                                 unlist(labeledPeak)));
  labelsMatrix = labelsMatrix[order(labelsMatrix[,3],labelsMatrix[,4]),]; # order into isotopologue groups

  ### DATA CLEAN-UP: ###
  # Part I:
  # Remove duplicate pairs and pairs that are both labeled peaks for other base peaks
  numPutativeLabels = dim(labelsMatrix)[1];
  basePeaks = unique(labelsMatrix[,1]);
  numLabeledPeaks = length(basePeaks);
  outtakes = list();
  for ( i in 1:numPutativeLabels ) { 
    B = labelsMatrix[,2] == labelsMatrix[i,1];
    A = labelsMatrix[B,1];
    C = which(labelsMatrix[,1] %in% A);
    if ( any(labelsMatrix[C,2] == labelsMatrix[i,2]) ) {
      outtakes = c(outtakes,i);
      next;
    }
    if ( i < numPutativeLabels ) {
      A = (i+1):numPutativeLabels;
      if ( any(labelsMatrix[A,1] == labelsMatrix[i,1] && labelsMatrix[A,2] == labelsMatrix[i,2])) {
        outtakes = c(outtakes,i);
      }
    }
  }
  outtakes = unlist(outtakes);
  labelsMatrix = labelsMatrix[-outtakes,];

  # Part II:
  # Check for significant difference between labeling patterns in unlabeled and labeled samples;
  # record only those groups that are different at p-value = alpha cutoff
  numPutativeLabels = dim(labelsMatrix)[1];
  basePeaks = unique(labelsMatrix[,1]);
  numLabeledPeaks = length(basePeaks);

  # output lists:
  base = list(); # base peak associated with isotopologue group
  mz = list(); # m/z 
  ID = list(); # peak ID number in xcmsSet@groups; if user wants to plot EICs of isotopologues, input this to XCMS's getEIC routine
  RT = list(); # retention time
  absInt1 = list(); # mean absolute intensity of peak in unlabeled samples
  absInt2 = list(); # ibid in labeled samples
  relInt1 = list(); # mean relative intensity of peak in unlabeled samples
  relInt2 = list(); # ibid in labeled samples
  totInt1 = list(); # total absolute intensity of all peaks in isotopologue group in unlabeled samples
  totInt2 = list(); #ibid in labeled samples
  CVabsInt1 = list(); # coefficient of variation of total ion intensity of isotopologue group in unlabeled samples
  CVabsInt2 = list(); # ibid in labeled samples
  SDrelInt1 = list(); # std dev of relative intensity of each isotopologue peak in unlabeled samples
  SDrelInt2 = list(); # ibid in labeled samples  
  foldEnrichment = list(); # (relative intensity in labeled samples) / (relative intensity in unlabeled samples)
  pvalues = list(); # from t-test for difference between unlabeled and labeled samples 
  sampleIntensities = list(); # ion counts (of type "intChoice") for each peak in isotopologue group in every sample
  
  # process each isotopologue group:
  j = 1;
  
  for (i in 1:numLabeledPeaks) {
    a = basePeaks[i]; # groupID of base peak
    baseIntensities = c(intensities1[a,], intensities2[a,]);
    isotopologues = list();
    IDs = list();
    RTs = list();
    numisotopologues = 0;
    k = j;
    # count the number of labeled peaks in each group
    while (k <= numPutativeLabels && labelsMatrix[k,1] == a) {
      isotopologues = c(isotopologues, groupMzs[labelsMatrix[k,2]]);
      IDs = c(IDs, groupIDs[labelsMatrix[k,2]]);
      RTs = c(RTs, groupRTs[labelsMatrix[k,2]]);
      numisotopologues = numisotopologues + 1;
      k = k+1;
    }
    isotopologues = unlist(isotopologues);
    IDs= unlist(IDs);
    RTs = unlist(RTs); 

    # discard peaks of low intensities in unlabeled samples as candidate base peaks
    if ( mean(intensities1[a,]) < noiseCutoff) { 
      j = k;
      next;
    }
    
    # create and fill a matrix to contain intensities of labeled peaks 
    abs1 = list();
    abs2 = list();
    
    labeledIntensities = matrix(rep(0,numisotopologues*numSamples),nrow = numisotopologues, ncol = numSamples);
    
    for (l in 1:numisotopologues) {
      b = labelsMatrix[j+l-1,2]; # groupID of labeled peak
      labeledIntensities[l,] = cbind(intensities1[b,,drop=FALSE], intensities2[b,,drop=FALSE]);
      abs1 = c(abs1, mean(intensities1[b,]));
      abs2 = c(abs2, mean(intensities2[b,]));
    }
    abs1 = unlist(abs1);
    abs2 = unlist(abs2);

    # remove redundantly called isotopologues by keeping only those with the smallest ppm error from expected m/z
    if (numisotopologues != length(unique(round(isotopologues)))) {
      M0 = round(groupMzs[a]);
      isos = round(isotopologues);
      reduced = unique(isos);
      numUniqIsos = length(reduced);
      outtakes = list();
      for (r in 1:numUniqIsos) {
        q = which(isos == reduced[r]);
        if (length(q) > 1 ) {
          massdefect = iMD * (reduced[r] - M0);
          delta = abs(groupMzs[IDs[q]] - groupMzs[a] - massdefect);
          outtakes = c(outtakes, q[which(delta != min(delta))]);          
        }
      }
      if (length(outtakes) >0) {
        outtakes = unlist(outtakes);
        isotopologues = isotopologues[-outtakes];
        IDs = IDs[-outtakes];
        RTs = RTs[-outtakes];
        numisotopologues = length(isotopologues);
        abs1 = abs1[-outtakes];
        abs2 = abs2[-outtakes];
        labeledIntensities = labeledIntensities[-outtakes,,drop=FALSE];
      }
    }

    # remove isotopologues whose intensities do not decrease monotonically in unlabeled samples from M0 to Mx   
    if (!compareOnlyDistros && monotonicityTol) {
      meanMprevUL = mean(intensities1[a,]); # mean intensity of M0 in unlabeled samples
      outtakes = list();
      for (l in 1:numisotopologues) {
        if (l ==1 && abs1[l] > (1+monotonicityTol)*meanMprevUL) {
          outtakes = c(outtakes,l);
        }
        else if (l > 1 && abs1[l] > (1+monotonicityTol)*meanMprevUL && round(isotopologues[l] - isotopologues[l-1]) > 1) {
          outtakes = c(outtakes,l);
        }
        else {
          meanMprevUL = abs1[l];
        }
      }
      outtakes = unlist(outtakes);
      if (length(outtakes) > 0) {
        abs1 = abs1[-outtakes];
        abs2 = abs2[-outtakes];
        labeledIntensities = labeledIntensities[-outtakes,,drop=FALSE];
        isotopologues = isotopologues[-outtakes];
        IDs= IDs[-outtakes];
        RTs = RTs[-outtakes];
      }
    }
    
    isotopologues = c(groupMzs[a], isotopologues);
    IDs= c(groupIDs[a], IDs);
    RTs = c(groupRTs[a], RTs); 
    allIntensities = rbind(baseIntensities, labeledIntensities);
    abs1 = c(mean(intensities1[a,]), abs1);
    abs2 = c(mean(intensities2[a,]), abs2);

    numisotopologues = length(isotopologues);
    sumIntensities = colSums(allIntensities);
    tot1 = mean(sumIntensities[1:dim(intensities1)[2]]); 
    tot2 = mean(sumIntensities[(dim(intensities1)[2]+1):numSamples]);
    cv1 = sd(sumIntensities[1:dim(intensities1)[2]])/tot1;
    cv2 = sd(sumIntensities[(dim(intensities1)[2]+1):numSamples])/tot2;

 
    # obtain relative intensities
    groupIntensities = allIntensities/ 
      matrix(rep(sumIntensities,numisotopologues),nrow = numisotopologues, byrow = TRUE);
    
    gI1 = groupIntensities[,1:dim(intensities1)[2], drop = FALSE];
    gI2 = groupIntensities[,(dim(intensities1)[2]+1):numSamples, drop = FALSE];
    
    # discard any samples with no signal for the isotopologue group    
    gI1 = gI1[,colSums(is.na(gI1))==0, drop = FALSE];
    gI2 = gI2[,colSums(is.na(gI2))==0, drop = FALSE];

    # discard whole group if majority of samples lack signal for the group
    if (dim(gI1)[2] < dim(intensities1)[2]/2 || dim(gI2)[2] < dim(intensities2)[2]/2) {
      j = k;
      next;
    }
    
    # calculate mean relative intensities
    rel1 = rowMeans(gI1);
    rel2 = rowMeans(gI2);
    sd1 = apply(gI1, 1, sd);
    sd2 = apply(gI2, 1, sd);

    enrichRatios = rel2/rel1;
    
    # discard whole group if enrichment is higher in unlabeled samples than in labeled ones, by a factor of 1+enrichTol
    if (!compareOnlyDistros) {
      if (enrichRatios[1] > (1+enrichTol)) {
        j = k;
        next;
      }
    }
    
    if (!singleSample) {
      pvalue = list();
      for (l in 1:numisotopologues) {
        if( all(gI1[l,] == 1) && all(gI2[l,] == 0) || is.infinite(enrichRatios[l]) ) {
          pvalue = c(pvalue,0);
        }
        else {
          T = try(t.test(gI1[l,], gI2[l,], var.equal = varEq), silent = TRUE);
          if (class(T) == "try-error") {
            pvalue = c(pvalue, 1);
            break;
          }
          else {
            pvalue = c(pvalue, T$p.value);
          }
        }
      }     

      # store data if at least one labeled isotopologue is enriched in labeled samples
      if ( any(unlist(pvalue) < alpha) && !any(unlist(pvalue) == 1) ) {
        base = c(base, groupMzs[a]);
        mz = c(mz, list(isotopologues));
        ID = c(ID, list(IDs)); 
        RT = c(RT, list(RTs)); 
        absInt1 = c(absInt1, list(abs1));
        absInt2 = c(absInt2, list(abs2)); 
        relInt1 = c(relInt1, list(rel1)); 
        relInt2 = c(relInt2, list(rel2));
        CVabsInt1 = c(CVabsInt1, cv1);
        CVabsInt2 = c(CVabsInt2, cv2); 
        totInt1 = c(totInt1, tot1);
        totInt2 = c(totInt2, tot2);
        SDrelInt1 = c(SDrelInt1, list(sd1));
        SDrelInt2 = c(SDrelInt2, list(sd2)); 
        foldEnrichment = c(foldEnrichment, list(enrichRatios)); 
        pvalues = c(pvalues, list(unlist(pvalue)));
        sampleIntensities = c(sampleIntensities,list(allIntensities));
      }
    }
    else { # if only single replicate is available for unlabeled and labeled samples, record all isotopologue groups
      deltaSpec = sum(abs(rel1[1:numisotopologues-1] - rel2[1:numisotopologues-1])); # differences in relative intensity profiles
      base = c(base, groupMzs[a]);
      mz = c(mz, list(isotopologues));
      ID = c(ID, list(IDs)); 
      RT = c(RT, list(RTs)); 
      absInt1 = c(absInt1, list(abs1));
      absInt2 = c(absInt2, list(abs2)); 
      relInt1 = c(relInt1, list(rel1)); 
      relInt2 = c(relInt2, list(rel2));
      CVabsInt1 = c(CVabsInt1, cv1);
      CVabsInt2 = c(CVabsInt2, cv2);
      totInt1 = c(totInt1, tot1);
      totInt2 = c(totInt2, tot2);
      SDrelInt1 = c(SDrelInt1, list(sd1));
      SDrelInt2 = c(SDrelInt2, list(sd2));
      foldEnrichment = c(foldEnrichment, list(enrichRatios)); 
      pvalues = c(pvalues, deltaSpec); # not real p-values
      sampleIntensities = c(sampleIntensities,list(allIntensities));
    }
    j = k;
  }

  labelsData = list(compound = base, isotopologue = mz, groupID = ID, rt = RT, meanAbsU = absInt1, totalAbsU = totInt1, cvTotalU = CVabsInt1,
                    meanAbsL = absInt2, totalAbsL = totInt2, cvTotalL = CVabsInt2, meanRelU = relInt1, meanRelL = relInt2, p_value = pvalues, enrichmentLvsU = foldEnrichment,
                    sdRelU = SDrelInt1, sdRelL = SDrelInt2, sampleData = sampleIntensities);
  
  return(labelsData);  
}

# outputs comparison of labeling reports of samples from conditions 1 and 2
getIsoDiffReport <- function(labelsData1, labelsData2, condition1, condition2, classes1, classes2, labeledSamples, varEq = FALSE, singleSample = FALSE) {
  # Variables:
  # labelsData1, *2: outputs of getIsoLabelReport() for samples from conditions 1 and 2
  # classes1, *2: character vectors designating whether each sample in labelsDataX is unlabeled or labeled
  #               (e.g. c("C12", "C12", "C13", "C13") if first two samples are unlabeled and second two are labeled)
  # labeledSamples: character variable name of labeled samples (e.g. "C13")
  # The variable names in classes1, classes2, and labeledSamples do not have to match the phenoData variable in an XCMS
  # object; they only have to match each other.
  # condition1, *2: character variables naming conditions 1 and 2 (e.g. "WT" and "KO")
  # singleSample: Boolean indicating whether only single samples were used to generate getIsoLabelReport for both sample classes
  
  # unpack the labeling reports:
  compounds1 = unlist(labelsData1$compound);
  compounds2 = unlist(labelsData2$compound);
  
  rt1 = labelsData1$rt;
  rt2 = labelsData2$rt;
   
#   enrich1 = labelsData1$enrichmentLvsU;
#   enrich2 = labelsData2$enrichmentLvsU;
  
  rel1U = labelsData1$meanRelU;
  rel1L = labelsData1$meanRelL;

#  sdrel1U = labelsData1$sdRelU;
  sdrel1L = labelsData1$sdRelL;
  
  rel2U = labelsData2$meanRelU;
  rel2L = labelsData2$meanRelL;
  
#  sdrel2U = labelsData2$sdRelU;
  sdrel2L = labelsData2$sdRelL;
  
#  abs1U = labelsData1$meanAbsU;
  abs1L = labelsData1$meanAbsL;
  
#  abs2U = labelsData2$meanAbsU;
  abs2L = labelsData2$meanAbsL;
    
  pattern1 = labelsData1$isotopologue;
  pattern2 = labelsData2$isotopologue;
  
  sampData1 = labelsData1$sampleData;
  sampData2 = labelsData2$sampleData;
  
  n1 = length(compounds1);
  n2 = length(compounds2);
  n2hasMatched = rep(NA,n2);
  n1_labeled = length(which(classes1 == labeledSamples));
  n2_labeled = length(which(classes2 == labeledSamples));

  # output lists:
  base = list();
  class1 = list();
  class2 = list();
  RT1 = list();
  RT2 = list();
#  enrichment1 = list();
#  enrichment2 = list();
  relints1u = list();
  relints1l = list();
  relints2u = list();
  relints2l = list();
#  SDrelints1u = list();
  SDrelints1l = list();
#  Sdrelints2u = list();
  SDrelints2l = list();

#  absints1u = list();
  absints1l = list();
#  absints2u = list();
  absints2l = list();

  pvalues = list(); 

  # search labeling report 1 against report 2:
  for (i in 1:n1) {
    match1found = 0;
    searchCompounds2 = compounds2[is.na(n2hasMatched)];
    len2 = length(searchCompounds2);
    for (j in 1:len2) {
      if ( compounds1[i] == searchCompounds2[j] ) { # if an isotopologue group is found in both conditions
        # if isotopologue list for this particular compound is not the same between the two groups:
        if ( length(pattern1[[i]]) != length(pattern2[is.na(n2hasMatched)][[j]])
             || any(pattern1[[i]] != pattern2[is.na(n2hasMatched)][[j]]) ) {
          p1 = pattern1[[i]];
          p2 = pattern2[is.na(n2hasMatched)][[j]];
          consensus = unique(c(p1,p2));
          consensus = consensus[order(consensus)];
          c1 = list();
          c2 = list();
          Ints1 = list();
          Ints2 = list();
#          e1 = list();
#          e2 = list();
          rT1 = list();
          rT2 = list();
          RI1u = list();
          RI1l = list();
          RI2u = list();
          RI2l = list();
          sdr1l = list();
          sdr2l = list();
#          AI1u = list();
          AI1l = list();
#          AI2u = list();
          AI2l = list();
          for ( k in 1:length(consensus) ) {
            P1 = match(consensus[k],p1);
            P2 = match(consensus[k],p2);
            if ( !is.na(P1) ) {
              Ints1 = c(Ints1,list(sampData1[[i]][P1,]));
              c1 = c(c1, pattern1[[i]][P1]);
 #             e1 = c(e1, enrich1[[i]][P1]);
              rT1 = c(rT1, rt1[[i]][P1]);
              RI1u = c(RI1u, rel1U[[i]][P1]);
              RI1l = c(RI1l, rel1L[[i]][P1]);
              sdr1l = c(sdr1l, sdrel1L[[i]][P1]);
#              AI1u = c(AI1u, abs1U[[i]][P1]);
              AI1l = c(AI1l, abs1L[[i]][P1]);
            }
            else {
              Ints1 = c(Ints1, list(rep(0, dim(sampData1[[i]])[2])));
              c1 = c(c1, NA);
#              e1 = c(e1, NA);
              rT1 = c(rT1, NA);
              RI1u = c(RI1u, NA);
              RI1l = c(RI1l, NA);
              sdr1l = c(sdr1l, NA);
#              AI1u = c(AI1u, NA);
              AI1l = c(AI1l, NA);
            }
            if ( !is.na(P2) ) {
              Ints2 = c(Ints2,list(sampData2[is.na(n2hasMatched)][[j]][P2,]));
              c2 = c(c2, pattern2[is.na(n2hasMatched)][[j]][P2]);
#              e2 = c(e2, enrich2[is.na(n2hasMatched)][[j]][P2]);
              rT2 = c(rT2, rt2[is.na(n2hasMatched)][[j]][P2]);
              RI2u = c(RI2u, rel2U[is.na(n2hasMatched)][[j]][P2]);
              RI2l = c(RI2l, rel2L[is.na(n2hasMatched)][[j]][P2]);
              sdr2l = c(sdr2l, sdrel2L[is.na(n2hasMatched)][[j]][P2]);
#              AI2u = c(AI2u, abs2U[is.na(n2hasMatched)][[j]][P2]);
              AI2l = c(AI2l, abs2L[is.na(n2hasMatched)][[j]][P2]);
            }
            else {
              Ints2 = c(Ints2, list(rep(0, dim(sampData2[is.na(n2hasMatched)][[j]])[2])));
              c2 = c(c2, NA);
#              e2 = c(e2, NA);
              rT2 = c(rT2, NA);
              RI2u = c(RI2u, NA);
              RI2l = c(RI2l, NA);
              sdr2l = c(sdr2l, NA);
#              AI2u = c(AI2u, NA);
              AI2l = c(AI2l, NA);
            }
          }
          
          class1 = c(class1, list(unlist(c1)));
          class2 = c(class2, list(unlist(c2)));
          Ints1 = matrix(unlist(Ints1), nrow = length(consensus), ncol = dim(sampData1[[i]])[2], byrow = TRUE);
          Ints2 = matrix(unlist(Ints2), nrow = length(consensus), ncol = dim(sampData2[is.na(n2hasMatched)][[j]])[2], byrow = TRUE);
#          enrichment1 = c(enrichment1, list(unlist(e1)));
#          enrichment2 = c(enrichment2, list(unlist(e2)));
          RT1 = c(RT1, list(unlist(rT1)));
          RT2 = c(RT2, list(unlist(rT2)));
          relints1u = c(relints1u, list(unlist(RI1u)));
          relints1l = c(relints1l, list(unlist(RI1l)));
          relints2u = c(relints2u, list(unlist(RI2u)));
          relints2l = c(relints2l, list(unlist(RI2l)));
          SDrelints1l = c(SDrelints1l, list(unlist(sdr1l)));
          SDrelints2l = c(SDrelints2l, list(unlist(sdr2l)));
#          absints1u = c(absints1u, list(unlist(AI1u)));
          absints1l = c(absints1l, list(unlist(AI1l)));
#          absints2u = c(absints2u, list(unlist(AI2u)));
          absints2l = c(absints2l, list(unlist(AI2l)));
        }
        else {
          class1 = c(class1, pattern1[i]);
          class2 = c(class2, pattern2[is.na(n2hasMatched)][j]);
          Ints1 = sampData1[[i]];
          Ints2 = sampData2[is.na(n2hasMatched)][[j]];
#          enrichment1 = c(enrichment1, enrich1[i]);
#          enrichment2 = c(enrichment2, enrich2[is.na(n2hasMatched)][j]);
          RT1 = c(RT1, rt1[i]);
          RT2 = c(RT2, rt2[is.na(n2hasMatched)][j]);
          relints1u = c(relints1u, rel1U[i]);
          relints1l = c(relints1l, rel1L[i]);
          relints2u = c(relints2u, rel2U[is.na(n2hasMatched)][j]);
          relints2l = c(relints2l, rel2L[is.na(n2hasMatched)][j]);
          SDrelints1l = c(SDrelints1l, sdrel1L[i]);
          SDrelints2l = c(SDrelints2l, sdrel2L[is.na(n2hasMatched)][j]);          
#          absints1u = c(absints1u, abs1U[i]);
          absints1l = c(absints1l, abs1L[i]);
#          absints2u = c(absints2u, abs2U[is.na(n2hasMatched)][j]);
          absints2l = c(absints2l, abs2L[is.na(n2hasMatched)][j]);        
        }
        
        # get relative intensities

        absInts1_labeled = Ints1[,which(classes1==labeledSamples), drop = FALSE];
        absInts2_labeled = Ints2[,which(classes2==labeledSamples), drop = FALSE];
        
        relInts1_labeled = absInts1_labeled / 
          matrix(rep(colSums(absInts1_labeled),dim(Ints1)[1]),nrow = dim(Ints1)[1], byrow = TRUE);
        relInts2_labeled = absInts2_labeled / 
          matrix(rep(colSums(absInts2_labeled),dim(Ints2)[1]),nrow = dim(Ints2)[1], byrow = TRUE);
        
        # calculate p-values associated with the differences between conditions 1 and 2 in
        # relative intensity of each isotopologue within group:
        if (singleSample == FALSE) {
          pvalue = list();
          for (l in 1:dim(relInts1_labeled)[1]) {
            T = try(t.test(relInts1_labeled[l,], relInts2_labeled[l,], var.equal = varEq), silent = TRUE);
            if (class(T) == "try-error" || is.na(T$p.value) ) {
              pvalue = c(pvalue, 1);
            }
            else {
              pvalue = c(pvalue, T$p.value);
            }
          }
          pvalues = c(pvalues, list(unlist(pvalue)));
        }
        else {
          RI1l = unlist(relints1l[[length(relints1l)]]);
          RI2l = unlist(relints2l[[length(relints2l)]]);
          RI1l[which(is.na(RI1l))] = 0;
          RI2l[which(is.na(RI2l))] = 0; 
          pvalue = sum(abs(RI1l - RI2l));
          pvalues = c(pvalues, pvalue);
        }

        
        base = c(base, compounds1[i]);
        n2hasMatched[is.na(n2hasMatched)][j] = 1;
        match1found = 1;
        break;
      }
    }
    # if isotopologue group only appears in condition 1
    if (!match1found) {
      base = c(base, compounds1[i]);
      class1 = c(class1, pattern1[i]);
      class2 = c(class2, -1);
      RT1 = c(RT1, rt1[i]);
      RT2 = c(RT2, -1);
#      enrichment1 = c(enrichment1, enrich1[i]);
#      enrichment2 = c(enrichment2, -1);
      pvalues = c(pvalues, 0);
      relints1u = c(relints1u, rel1U[i]);
      relints1l = c(relints1l, rel1L[i]);
      relints2u = c(relints2u, -1);
      relints2l = c(relints2l, -1);
      SDrelints1l = c(SDrelints1l, sdrel1L[i]);
      SDrelints2l = c(SDrelints2l, -1);
#      absints1u = c(absints1u, abs1U[i]);
      absints1l = c(absints1l, abs1L[i]);
#      absints2u = c(absints2u, -1);
      absints2l = c(absints2l, -1);  
    }
  }
 
  # process isotopologue groups that only appear in condition 2
  unmatchedCompounds2 = compounds2[is.na(n2hasMatched)];
  len2 = length(unmatchedCompounds2);
  for (i in 1:len2) {
    base = c(base, unmatchedCompounds2[i]);
    class1 = c(class1, -1);
    class2 = c(class2, pattern2[is.na(n2hasMatched)][i]);
    RT1 = c(RT1, -1);
    RT2 = c(RT2, rt2[is.na(n2hasMatched)][i]);
#    enrichment1 = c(enrichment1, -1);
#    enrichment2 = c(enrichment2, enrich2[is.na(n2hasMatched)][i]);
    pvalues = c(pvalues, 0);
    relints1u = c(relints1u, -1);
    relints1l = c(relints1l, -1);
    relints2u = c(relints2u, rel2U[is.na(n2hasMatched)][i]);
    relints2l = c(relints2l, rel2L[is.na(n2hasMatched)][i]);
    SDrelints1l = c(SDrelints1l, -1);
    SDrelints2l = c(SDrelints2l, sdrel2L[is.na(n2hasMatched)][i]);
#    absints1u = c(absints1u, -1);
    absints1l = c(absints1l, -1);
#    absints2u = c(absints2u, abs2U[is.na(n2hasMatched)][i]);
    absints2l = c(absints2l, abs2L[is.na(n2hasMatched)][i]);
  }

  isoDiffReport = list(class1[order(unlist(base))], class2[order(unlist(base))], 
                       RT1[order(unlist(base))], RT2[order(unlist(base))], 
#                       enrichment1[order(unlist(base))], enrichment2[order(unlist(base))], 
                       relints1u[order(unlist(base))], relints1l[order(unlist(base))], 
                       relints2u[order(unlist(base))], relints2l[order(unlist(base))], 
#                       absints1u[order(unlist(base))], absints1l[order(unlist(base))],
#                       absints2u[order(unlist(base))], absints2l[order(unlist(base))], 
                       pvalues[order(unlist(base))],
                       SDrelints1l[order(unlist(base))], SDrelints2l[order(unlist(base))],
                       absints1l[order(unlist(base))], absints2l[order(unlist(base))]);
#  names(isoDiffReport) = c(condition1, condition2, "rt1", "rt2", "isotopeEnrichment1", "isotopeEnrichment2", "relInts1U", "relInts1L", "relInts2U", "relInts2L",
#                           "absInts1U", "absInts1L", "absInts2U", "absInts2L","p_value");
  names(isoDiffReport) = c(condition1, condition2, "rt1", "rt2", "relInts1U", "relInts1L", "relInts2U", "relInts2L",
                         "p_value", "sdrelInts1L", "sdrelInts2L", "absInts1L", "absInts2L");  
  return(isoDiffReport); 
}

# prints either getIsoLabelReport() or getIsoDiffReport() outputs to tab-delimited text file named "outputfile"
printIsoListOutputs <- function(listReport, outputfile) {
  colNames = names(listReport);
  nblocks = length(colNames);
  colNames = c("", colNames);
  write(colNames, file = outputfile, ncolumns = nblocks+1, append = FALSE, sep = "\t");
  
  nrows = length(listReport[[1]]);
  for (i in 1:nrows) {
    maxRows = 1;
    ncols = 0;
    for (j in 1:nblocks) {
      cell = listReport[[j]][[i]];
      if ( length(cell) > maxRows && class(cell) == "numeric") {
        maxRows = length(cell);
        ncols = ncols + 1;
      }
      else if ( class(cell) == "matrix" ) {
        ncols = ncols + dim(cell)[2];
      }
      else {
        ncols = ncols+1;
      }
    }
    rowMatrix = matrix(rep(NA, (maxRows+1)*ncols),nrow = maxRows+1);
    for (j in 1:maxRows) {
      ncol = 1;
      for (k in 1:nblocks) {
        cell = listReport[[k]][[i]]; 
        if ( length(cell) == 1 && j == 1 && cell != -1 ) {
          rowMatrix[j, ncol] = cell;
          ncol = ncol + 1;
        }
        else if ( class(cell) == "numeric" || class(cell) == "integer" || class(cell) == "character") {
          if ( !is.na(cell[j]) && cell[j] != -1 ) {
            rowMatrix[j, ncol] = cell[j];
          }
          ncol = ncol + 1;
        }
        else if ( class(cell) == "matrix" ) {
          rowMatrix[j, ncol:(ncol+dim(cell)[2]-1)] = cell[j,];
          ncol = ncol+dim(cell)[2];
        }
      }
    }
    nr = dim(rowMatrix)[1];
    rowMatrix = cbind(c(i,rep(NA,nr-1)),rowMatrix);
    write.table(rowMatrix,outputfile,na = "",row.names = FALSE, col.names = FALSE, 
                sep = "\t", append = TRUE);
  }
}

# plots isotopologue distributions from a label report to PDF named outputfile; 12 plots/page
plotLabelReport <- function(isoLabelReport, intOption = "rel", classes, labeledSamples, outputfile) {
  pdf(outputfile, width = 8.5, height = 11);
  numGroups = length(isoLabelReport[[1]]);
  par(mfrow=c(4,3));
  for (i in 1:numGroups) {
    if (i %% 12 == 1) {
      par(mfrow=c(4,3));
    }
    if (intOption == "rel") {
      toPlot = cbind(isoLabelReport$meanRelU[[i]], isoLabelReport$meanRelL[[i]]);
      sds = cbind(isoLabelReport$sdRelU[[i]], isoLabelReport$sdRelL[[i]]);
      YLIM = c(0,1);
    }
    else if (intOption == "abs") {
      toPlot = cbind(isoLabelReport$meanAbsU[[i]], isoLabelReport$meanAbsL[[i]]);
      sds = cbind(isoLabelReport$cvTotalU[[i]]*isoLabelReport$meanAbsU[[i]], isoLabelReport$cvTotalL[[i]]*isoLabelReport$meanAbsL[[i]]);
      YLIM = NULL;
      totalPoolsU = colSums(isoLabelReport$sampleData[[i]][,which(classes != labeledSamples)]);
      totalPoolsL = colSums(isoLabelReport$sampleData[[i]][,which(classes == labeledSamples)]);
      T = t.test(totalPoolsU, totalPoolsL);
      pval = round(T$p.value, 3);
    }
    else {
      print("Error. Specify intOption as rel or abs.");
      break;
    }
    numPeaks = dim(toPlot)[1];
    
    rownames(toPlot) = as.character(round(isoLabelReport$isotopologue[[i]], 4));
    colnames(toPlot) = c("U", "L");
      
    if (intOption == "rel") {
      mp = barplot(t(toPlot), beside = TRUE, legend.text = c("U", "L"), 
              main = paste("m/z = ", rownames(toPlot)[1], "\n", "RT = ", isoLabelReport$rt[[i]][1]), 
              axisnames = FALSE, ylim = YLIM);
      text(seq(1.5, (numPeaks-1)*3+1.5, by = 3), par("usr")[3]-0.2, xpd = TRUE, labels = rownames(toPlot), srt = 45);
      
      segments(mp, t(toPlot) - t(sds), mp, t(toPlot) + t(sds), lwd = 2 );
    }
    else {
      barplot(toPlot, legend.text = rownames(toPlot), 
              main = paste("m/z = ", rownames(toPlot)[1], "\n", "RT = ", isoLabelReport$rt[[i]][1], "\n", "p = ", pval));
    }
  }
  dev.off();
}

# plots isotopologue distributions (represented as relative intensities) from isoDiff report to PDF named outputfile; 12 plots/page
plotIsoDiffReport <- function(isoDiffReport, xcmsSet, intChoice, sampleNames1, sampleNames2, 
                              labelReport1, labelReport2, classes1, classes2, labeledSamples, 
                              isotopeMassDifference, outputfile) {
  pdf(outputfile, width = 8.5, height = 11);
  numGroups = length(isoDiffReport[[1]]);
  par(mfrow=c(4,3));
  ints = groupval(xcmsSet, "medret", intChoice);
  ints1 = ints[,which(colnames(ints) %in% sampleNames1[which(classes1 == labeledSamples)])];
  ints2 = ints[,which(colnames(ints) %in% sampleNames2[which(classes2 == labeledSamples)])];
  ints1[is.na(ints1)] = 0;
  ints2[is.na(ints2)] = 0;
  for (i in 1:numGroups) {
    if (i %% 12 == 1) {
      par(mfrow=c(4,3));
    }
    if (isoDiffReport$relInts1L[[i]][1] != -1 &&
          isoDiffReport$relInts2L[[i]][1] != -1) {
      toPlot = cbind(isoDiffReport$relInts1L[[i]], isoDiffReport$relInts2L[[i]]);
      sds = cbind(isoDiffReport$sdrelInts1L[[i]], isoDiffReport$sdrelInts2L[[i]]);
      numPeaks = dim(toPlot)[1];
      a = isoDiffReport[[1]][[i]];
      b = isoDiffReport[[2]][[i]];
      if ( NA %in% a ) {
        mz = round(b[1], 4);
        isos = round((b-b[1])/isotopeMassDifference);
        rownames(toPlot) = as.character(isos);
      }
      else {
        mz = round(a[1], 4);
        isos = round((a-a[1])/isotopeMassDifference);
        rownames(toPlot) = as.character(isos);
      }
      mp = barplot(t(toPlot), beside = TRUE,legend.text = names(isoDiffReport)[1:2], 
              main = paste("m/z = ", mz, "\n", "RT = ", isoDiffReport$rt1[[i]][1]));
#       mp = barplot(t(toPlot), beside = TRUE, legend.text = names(isoDiffReport)[1:2], 
#               main = paste("m/z = ", rownames(toPlot)[1], "\n", "RT = ", isoDiffReport$rt1[[i]][1]), 
#               axisnames = FALSE, ylim = c(0,1));
#       text(seq(1.5, (numPeaks-1)*3+1.5, by = 3), par("usr")[3]-0.1, xpd = TRUE, labels = rownames(toPlot), srt = 45);
      sds[which(is.na(sds))] = 0;
      segments(mp, t(toPlot) - t(sds), mp, t(toPlot) + t(sds), lwd = 2 );      
    }
    else if (isoDiffReport$relInts1L[[i]][1] == -1) {
      a = length(isoDiffReport$relInts2L[[i]]);
      b = isoDiffReport[[2]][[i]];
      ind = which(unlist(labelReport2$compound) == b[1]); 
      abs = ints1[labelReport2$groupID[[ind]],];
      tots = colSums(abs);
      rel = abs/matrix(rep(tots,a), byrow = TRUE, nrow = a);
      meanRel = rowMeans(rel);
      sd = apply(rel, 1, sd);
      mz = round(isoDiffReport[[2]][[i]][1], 4);
      toPlot = cbind(meanRel, isoDiffReport$relInts2L[[i]]);
      sds = cbind(sd, isoDiffReport$sdrelInts2L[[i]]);
      numPeaks = dim(toPlot)[1];
      isos = round((isoDiffReport[[2]][[i]]-isoDiffReport[[2]][[i]][1])/isotopeMassDifference);
      rownames(toPlot) = as.character(isos);
      mp = barplot(t(toPlot), beside = TRUE, legend.text = names(isoDiffReport)[1:2],
              main = paste("m/z = ", mz, "\n", "RT = ", isoDiffReport$rt2[[i]][1]));
              #               main = paste("m/z = ", rownames(toPlot)[1], "\n", "RT = ", isoDiffReport$rt2[[i]][1]), 
#       mp = barplot(t(toPlot), beside = TRUE, legend.text = names(isoDiffReport)[1:2],
#               main = paste("m/z = ", rownames(toPlot)[1], "\n", "RT = ", isoDiffReport$rt2[[i]][1]), 
#               axisnames = FALSE, ylim = c(0,1));
#       text(seq(1.5, (numPeaks-1)*3+1.5, by = 3), par("usr")[3]-0.1, xpd = TRUE, labels = rownames(toPlot), srt = 45);
      sds[which(is.na(sds))] = 0;
      segments(mp, t(toPlot) - t(sds), mp, t(toPlot) + t(sds), lwd = 2 );   
    }
    else {
      a = length(isoDiffReport$relInts1L[[i]]);
      b = isoDiffReport[[1]][[i]];
      ind = which(unlist(labelReport1$compound) == b[1]); 
      abs = ints2[labelReport1$groupID[[ind]],];
      tots = colSums(abs);
      rel = abs/matrix(rep(tots,a), byrow = TRUE, nrow = a);
      meanRel = rowMeans(rel);
      sd = apply(rel, 1, sd);
      mz = round(isoDiffReport[[1]][[i]][1], 4);
      toPlot = cbind(isoDiffReport$relInts1L[[i]],meanRel);
      sds = cbind(isoDiffReport$sdrelInts1L[[i]], sd); 
      numPeaks = dim(toPlot)[1];
      isos = round((isoDiffReport[[1]][[i]]-isoDiffReport[[1]][[i]][1])/isotopeMassDifference);
      rownames(toPlot) = as.character(isos);
      mp = barplot(t(toPlot), beside = TRUE, legend.text = names(isoDiffReport)[1:2],
               main = paste("m/z = ", mz, "\n", "RT = ", isoDiffReport$rt1[[i]][1])); 
#       mp = barplot(t(toPlot), beside = TRUE, legend.text = names(isoDiffReport)[1:2],
#               main = paste("m/z = ", rownames(toPlot)[1], "\n", "RT = ", isoDiffReport$rt1[[i]][1]), 
#               axisnames = FALSE, ylim = c(0,1));
#       text(seq(1.5, (numPeaks-1)*3+1.5, by= 3), par("usr")[3]-0.1, xpd = TRUE, labels = rownames(toPlot), srt = 45);
      sds[which(is.na(sds))] = 0;
      segments(mp, t(toPlot) - t(sds), mp, t(toPlot) + t(sds), lwd = 2 );      
    }
  }
  dev.off();
}

# plots total absolute intensities of every isotopologue group to PDF named outputfile; 12 plots/page
plotTotalIsoPools <- function(isoDiffReport, xcmsSet, intChoice, sampleNames1, sampleNames2, 
                              labelReport1, labelReport2, classes1, classes2, labeledSamples, 
                              outputfile) {
  pdf(outputfile, width = 8.5, height = 11);
  numGroups = length(isoDiffReport[[1]]);
  par(mfrow=c(4,3));
  ints = groupval(xcmsSet, "medret", intChoice);
  ints1 = ints[,which(colnames(ints) %in% sampleNames1[which(classes1 == labeledSamples)])];
  ints2 = ints[,which(colnames(ints) %in% sampleNames2[which(classes2 == labeledSamples)])];
  ints1[is.na(ints1)] = 0;
  ints2[is.na(ints2)] = 0;
  for (i in 1:numGroups) {
    if (i %% 12 == 1) {
      par(mfrow=c(4,3));
    }

    if (isoDiffReport$absInts1L[[i]][1] != -1 &&
          isoDiffReport$absInts2L[[i]][1] != -1) {
      toPlot = cbind(isoDiffReport$absInts1L[[i]], isoDiffReport$absInts2L[[i]]);
      toPlot[is.na(toPlot)] = 0;
      a = isoDiffReport[[1]][[i]];
      b = isoDiffReport[[2]][[i]];

      ind1 = which(unlist(labelReport1$compound) == a[1]);
      ind2 = which(unlist(labelReport2$compound) == b[1]);
      totalPools1 = colSums(labelReport1$sampleData[[ind1]][,which(classes1 == labeledSamples)]);
      totalPools2 = colSums(labelReport2$sampleData[[ind2]][,which(classes2 == labeledSamples)]);
      T = t.test(totalPools1, totalPools2);
      pval = round(T$p.value, 3);
      
      if ( NA %in% a ) {
        rownames(toPlot) = as.character(round(b-b[1]));
      }
      else {
        rownames(toPlot) = as.character(round(a-a[1]));
      }
      colnames(toPlot) = c(unlist(attributes(isoDiffReport))[1],unlist(attributes(isoDiffReport))[2]);
      barplot(toPlot, legend.text = rownames(toPlot), 
              main = paste("m/z =", round(b[1],4), ", ", "RT =", isoDiffReport$rt1[[i]][1], "\n", "p =", pval));

    }
    else if (isoDiffReport$absInts1L[[i]][1] == -1) {
      b = isoDiffReport[[2]][[i]];
      ind = which(unlist(labelReport2$compound) == b[1]); 
      a = rowMeans(ints1[labelReport2$groupID[[ind]],]);
      totalPools1 = colSums(ints1[labelReport2$groupID[[ind]],]);
      totalPools2 = colSums(labelReport2$sampleData[[ind]][,which(classes2 == labeledSamples)]);
      T = t.test(totalPools1, totalPools2);
      pval = round(T$p.value, 3);
      
      toPlot = cbind(a, isoDiffReport$absInts2L[[i]]);

      rownames(toPlot) = as.character(round(isoDiffReport[[2]][[i]]-isoDiffReport[[2]][[i]][1]));
      colnames(toPlot) = c(unlist(attributes(isoDiffReport))[1],unlist(attributes(isoDiffReport))[2]);
      barplot(toPlot, legend.text = rownames(toPlot),
              main = paste("m/z =", round(isoDiffReport[[2]][[i]][1],4), ", ", "RT =", isoDiffReport$rt2[[i]][1], "\n",
              "p =", pval, "\n", "(", colnames(toPlot)[1], "not significantly labeled )"));

    }
    else {
      b = isoDiffReport[[1]][[i]];
      ind = which(unlist(labelReport1$compound) == b[1]);
      a = rowMeans(ints2[labelReport1$groupID[[ind]],]);
      totalPools2 = colSums(ints2[labelReport1$groupID[[ind]],]);
      totalPools1 = colSums(labelReport1$sampleData[[ind]][,which(classes1 == labeledSamples)]);
      T = t.test(totalPools1, totalPools2);
      pval = round(T$p.value, 3);

      toPlot = cbind(isoDiffReport$absInts1L[[i]],a);
      rownames(toPlot) = as.character(round(isoDiffReport[[1]][[i]]-isoDiffReport[[1]][[i]][1]));
      colnames(toPlot) = c(unlist(attributes(isoDiffReport))[1],unlist(attributes(isoDiffReport))[2]);
      barplot(toPlot, legend.text = rownames(toPlot),
              main = paste("m/z =", round(isoDiffReport[[1]][[i]][1],4), ", ", "RT =", isoDiffReport$rt1[[i]][1], "\n",
              "p =", pval, "\n", "(", colnames(toPlot)[2], "not significantly labeled )"));

    }
  }
  dev.off();
}

# generates conventional XCMS-like diff report comparing peaks between two sample classes, e.g. "control", and "treatment" from an xcmsSet object
# configured for XC13MS (i.e. one in which the sample classes have to be demarcated as unlabeled and labeled)
miniDiffReport <- function(xcmsSet, class1sampNames, class2sampNames, varEq = FALSE, intChoice) {
  peakIntensities = groupval(xcmsSet, method = "medret", value = intChoice);
  if (intChoice == "intb") {
    peakIntensities[is.na(peakIntensities)] = 0;
  }
  peakIntensities1 = peakIntensities[, match(class1sampNames, colnames(peakIntensities))];
  peakIntensities2 = peakIntensities[, match(class2sampNames, colnames(peakIntensities))];
  means1 = rowMeans(peakIntensities1);
  means2 = rowMeans(peakIntensities2);
  groups = data.frame(xcmsSet@groups);
  npeaks = dim(peakIntensities)[1];
  pvalue = rep(0, npeaks);
  foldChange = rep(0, npeaks);
  for (i in 1:npeaks) { # i is also group number (for grabbing group info)
    T = t.test(peakIntensities1[i,], peakIntensities2[i,], var.equal = varEq);
    pvalue[i] = T$p.value;
    if ( means2[i] > means1[i] ) {
      foldChange[i] = means2[i] / means1[i];
    }
    else {
      foldChange[i] = -means1[i] / means2[i];
    }
  }
  foldChange = foldChange[order(pvalue)];
  peakIntensities1 = peakIntensities1[order(pvalue),];
  peakIntensities2 = peakIntensities2[order(pvalue),];
  means1 = means1[order(pvalue)];
  means2 = means2[order(pvalue)];
  groups = groups[order(pvalue), c(1,4)];
  pvalue = pvalue[order(pvalue)];
  out = cbind(pvalue, foldChange, groups, means1, means2, peakIntensities1, peakIntensities2);
  return(out);
}

# filters isoDiff report to report back only isotopologue groups that are different between sample classes
filterIsoDiffReport <- function(isoDiffReport, alpha, whichPeak = 2) {
  n = length(isoDiffReport[[1]]);
  outtakes = list();
  for (i in 1:n) {
    # filter 1: remove groups with p-values above alpha
    if (any(unlist(isoDiffReport$p_value[i]) == 1)) {
      outtakes = c(outtakes,i);
      next;
    }
    # if whichPeak == 1, only consider intensity differences of base isotopologue when calling significance
    if (whichPeak == 1 && isoDiffReport$p_value[i][[1]] > alpha ) {
      outtakes = c(outtakes,i);
      next;
    }
    # if whichPeak == 2, consider intensity differences of all isotopologues and call group 
    # as signficant if any of them are different (akin to default reporting from getIsoDiffReport())
    if (whichPeak == 2 && !any(unlist(isoDiffReport$p_value[i]) < alpha) ) {
      outtakes = c(outtakes,i);
      next;
    }
  }
  outtakes = unlist(outtakes);
  if (length(outtakes) == 0) {
    out = isoDiffReport;
  }
  else {
    out = lapply(isoDiffReport, "[", -outtakes);
  }
  return(out);
}