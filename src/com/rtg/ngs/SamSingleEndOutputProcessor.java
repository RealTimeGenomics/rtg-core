/*
 * Copyright (c) 2014. Real Time Genomics Limited.
 *
 * Use of this source code is bound by the Real Time Genomics Limited Software Licence Agreement
 * for Academic Non-commercial Research Purposes only.
 *
 * If you did not receive a license accompanying this file, a copy must first be obtained by email
 * from support@realtimegenomics.com.  On downloading, using and/or continuing to use this source
 * code you accept the terms of that license agreement and any amendments to those terms that may
 * be made from time to time by Real Time Genomics Limited.
 */
package com.rtg.ngs;


import java.io.File;
import java.io.IOException;
import java.util.Collections;

import com.rtg.launcher.HashingRegion;
import com.rtg.util.intervals.ReferenceRegions;
import com.rtg.calibrate.Calibrator;
import com.rtg.calibrate.CovariateEnum;
import com.rtg.calibrate.Recalibrate;
import com.rtg.index.hash.ngs.OutputProcessor;
import com.rtg.ngs.blocking.MapQScoringReadBlocker;
import com.rtg.ngs.blocking.ReadBlocker;
import com.rtg.ngs.blocking.ReadBlockerSync;
import com.rtg.pairedend.ReadStatusListener;
import com.rtg.reader.NamesInterface;
import com.rtg.sam.BamIndexer;
import com.rtg.tabix.TabixIndexer;
import com.rtg.tabix.UnindexableDataException;
import com.rtg.util.StringUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.Timer;

/**
 * Writes SAM output files for single end datasets.
 *
 */
public class SamSingleEndOutputProcessor extends AbstractMapOutputProcessor {

  private static final int MATCHED = ReadStatusTracker.MATCHED_FIRST | ReadStatusTracker.MATCHED_SECOND;
  private static final int BLOCKED = ReadStatusTracker.BLOCKED_FIRST | ReadStatusTracker.BLOCKED_SECOND;

  private final boolean mOutputUnmapped;
  private UptoNStore mTopN;
  private final ReadBlocker mFreqBlocker;


  /**
   * Create a new {@link SamSingleEndOutputProcessor}
   * @param param {@link NgsParams}
   * @param stats map to put statistics into
   * @param outputUnmapped  true if unmapped should be output
   * @throws IOException If an IO error
   */
  public SamSingleEndOutputProcessor(NgsParams param, MapStatistics stats, boolean outputUnmapped) throws IOException {
    super(param, stats, false, /** irrelevant */ false);
    //super((int) param.build().reader().numberSequences(), stats);
    final long numSequences = param.buildFirstParams().numberSequences();
    mOutputUnmapped = outputUnmapped;
    mFreqBlocker = new ReadBlockerSync((int) numSequences, param.readFreqThreshold(), "single end freq count");
    mTopN = new UptoNStoreSync(new TopNImplementation((int) numSequences, (int) param.searchParams().reader().numberSequences(), mParams.outputParams().filter().topN(), param.searchParams().reader().maxLength(), param.buildFirstParams().maxLength()));
  }

  @Override
  public void process(final long templateId, final String frame, final int readId, final int tStart, final int score, final int scoreIndel) {
    if (mFreqBlocker.isBlocked(readId)) {
      mUnmappedTracker.addStatus(readId, BLOCKED);
      return;
    }
    mFreqBlocker.increment(readId);
    mTopN.process(templateId, frame.startsWith("R"), readId, tStart, scoreIndel);
    mUnmappedTracker.addStatus(readId, MATCHED);
  }

  @Override
  public void finish() throws IOException {
    sortRegions();
    Diagnostic.developerLog(mTopN.histogram());
    final FilterConcatIntermediateFiles alignmentIntFiles = writeAlignments();

    final FilterConcatIntermediateFiles unmappedIntFiles;
    if (mOutputUnmapped) {
      unmappedIntFiles = writeUnmapped(!mParams.outputParams().unify(), false, false);
    } else {
      unmappedIntFiles = null;
    }
    if (mParams.outputParams().unify()) {
      whizBangUnify(unmappedIntFiles, alignmentIntFiles);
//      unifyAlignmentOutput(alignmentFiles);
    }
    mUnmappedTracker.calculateStatistics(false, false);
    mReportMerger.blendReportData().write(new File(mParams.outputParams().directory(), MapReportData.MAP_REPORT_FILE_NAME));

  }
  @Override
  public synchronized OutputProcessor threadClone(HashingRegion region) throws IOException {
    super.threadClone(region);
    if (region != HashingRegion.NONE) {
      return new ClippedOutputProcessor(this, region);
    } else {
      return this;
    }
  }

  @Override
  public void threadFinish() {
  }

  private static MatchResult organizeResults(final int[] readIdStatus, final UptoNStore uptoN) {
    final MatchResult results = new MatchResult(readIdStatus.length);
    for (int i = 0; i < readIdStatus.length; ++i) {
      //System.err.println(i + " READ ID STATUS = " + mReadIdStatus[i]);
      uptoN.setResults(results, i);

    }
    results.sort();
    return results;
  }


  FilterConcatIntermediateFiles writeAlignments() throws IOException {
    Diagnostic.userLog("Extracting hits for reads");
    final Timer outputTimer = new Timer("AlignmentOutput");
    outputTimer.start();
    final MatchResult results = organizeResults(mUnmappedTracker.mReadIdStatus, mTopN);
    mTopN = null;
    Collections.sort(mRegions);
    final FilterConcatIntermediateFiles files = getNonMatedFilterConcatIntermediateFiles(results, false, mRegions.toArray(new HashingRegion[mRegions.size()]));
    outputTimer.stop();
    outputTimer.log();
    return files;
    //indexSamFile(mParams, outFile);
  }

  @Override
  protected FilterConcatIntermediateFiles filterConcatNonMated(MapQScoringReadBlocker blockerLeft, MapQScoringReadBlocker blockerRight, File[] tempFiles, SingleEndTopRandomImplementation.HitRecord[] hitsToKeep, NamesInterface templateNames, File outFile) throws IOException {
    return new SingleEndMulticoreFilterConcat(mParams, mUnmappedTracker, blockerLeft, mFreqBlocker, hitsToKeep, templateNames, mReportMerger).filterConcat(tempFiles, outFile, mSharedResources.getHeader(), mParams.outputParams());
  }

  static boolean canIndex(final NgsParams params, final File outFile) {
    if (params.searchParams().reader().maxLength() > TabixIndexer.MAXIMUM_REFERENCE_LENGTH) {
      Diagnostic.warning("Cannot produce TABIX index for: " + outFile + " as maximum reference sequence length is exceeded.");
      return false;
    }
    return true;
  }

  // No better place for this to live
  static File indexSamFile(final NgsParams params, final File outFile, boolean hasHeader, int numSequences) throws IOException {
    File indexFile = null;
    if (params.blockCompressed() && params.outputParams().outputIndex() && canIndex(params, outFile)) {
      final Timer indexTimer = new Timer("IndexSam-" + outFile.getName());
      indexTimer.start();
      try {
        if (params.outputParams().bam()) {
          indexFile = new File(outFile.getParentFile(), outFile.getName() + BamIndexer.BAM_INDEX_EXTENSION);
          BamIndexer.saveBamIndexNoHeader(outFile, indexFile, hasHeader, numSequences);
        } else {
          indexFile = new File(outFile.getParentFile(), outFile.getName() + TabixIndexer.TABIX_EXTENSION);
          new TabixIndexer(outFile, indexFile).saveSamIndex();
        }
      } catch (final UnindexableDataException e) {
        Diagnostic.warning("Cannot produce TABIX index for: " + outFile + ": " + e.getMessage());
        indexFile = null;
      }
      indexTimer.stop();
      indexTimer.log();
    }
    return indexFile;
  }

  // No better place for this to live either
  static File calibrateUnmappedFile(final NgsParams params, final File outFile) throws IOException {
    if (params.outputParams().calibrate()) {
      final ReferenceRegions regions = params.outputParams().calibrateRegions();
      final Calibrator cal = new Calibrator(CovariateEnum.getCovariates(CovariateEnum.DEFAULT_COVARIATES, null), regions);
      if (regions != null) {
        cal.setSequenceLengths(Calibrator.getSequenceLengthMap(params.searchParams().reader(), regions));
      }
      final File calFile = new File(outFile.getParentFile(), outFile.getName() + Recalibrate.EXTENSION);
      cal.writeToFile(calFile);
      return calFile;
    }
    return null;
  }
  @Override
  public void close() {
  }

  private static class SingleEndMulticoreFilterConcat extends AbstractMulticoreFilterConcat {

    final MapQScoringReadBlocker mAsBlocker;
    final ReadBlocker mFreqBlocker;
    final long mReadIdOffset;
    final SingleEndTopRandomImplementation.HitRecord[] mHitsToKeep;
    final NamesInterface mTemplateNames;
    final ReadStatusListener mListener;
    final MapReportData.Merger mReportMerger;

    SingleEndMulticoreFilterConcat(NgsParams params, ReadStatusListener listener, MapQScoringReadBlocker asBlocker, ReadBlocker freqBlocker, SingleEndTopRandomImplementation.HitRecord[] hitsToKeep, NamesInterface templateNames, MapReportData.Merger reportMerger) {
      super(params);
      mListener = listener;
      mAsBlocker = asBlocker;
      mFreqBlocker = freqBlocker;
      mReadIdOffset = Math.max(0, params.buildFirstParams().readerRestriction().getStart());
      mHitsToKeep = hitsToKeep;
      mTemplateNames = templateNames;
      mReportMerger = reportMerger;
    }

    @Override
    protected AbstractSamResultsFilter makeFilter() {
      final String readGroupId = mParams.outputParams().readGroup() != null ? mParams.outputParams().readGroup().getReadGroupId() : null;
      final SingleEndSamResultsFilter filter = new SingleEndSamResultsFilter(mAsBlocker, mFreqBlocker, mListener, mReadIdOffset, mParams.buildFirstParams().reader().copy(), readGroupId, mParams.legacyCigars());
      filter.setBamOutput(mParams.outputParams().bam());
      filter.setMapReportData(mReportMerger);
      if (mHitsToKeep != null) {
        filter.setHitsToKeep(mHitsToKeep);
        filter.setTemplateNames(mTemplateNames);
      }
      if (mParams.outputParams().outputReadNames()) {
        try {
          assert mParams.buildSecondParams() == null;
          filter.setReadNames(mParams.buildFirstParams().reader().names());
        } catch (final IOException e) {
          filter.setReadNames(null);
          Diagnostic.warning("Failed to retrieve read names, using read id instead.");
        }
      }
      return filter;
    }
  }

  @Override
  public String toString() {
    final StringBuilder str = new StringBuilder();
    str.append("SamSingleEndOutputProcessor; topn= ").append(mTopN).append(StringUtils.LS);
    final long num = mParams.buildFirstParams().reader().numberSequences();
    str.append("numsequences= ").append(num).append(StringUtils.LS);
    str.append("temp files gzipped= " + true);
    for (int i = 0 ; i < num; ++i) {
      str.append("read= ").append(i).append(" stats= ").append(mUnmappedTracker.getXCAttribute(i, true)).append(StringUtils.LS);
    }
    return str.toString();
  }
}

