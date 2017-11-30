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
package com.rtg.protein;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.util.ArrayList;

import com.rtg.alignment.ActionsHelper;
import com.rtg.alignment.BidirectionalEditDistance;
import com.rtg.alignment.EditDistanceFactory;
import com.rtg.index.hash.ngs.OutputProcessor;
import com.rtg.mode.Frame;
import com.rtg.mode.ProteinScoringMatrix;
import com.rtg.mode.TranslatedFrame;
import com.rtg.ngs.MapStatistics;
import com.rtg.ngs.NgsParams;
import com.rtg.reader.SequencesReader;
import com.rtg.util.Environment;
import com.rtg.util.IntegerOrPercentage;
import com.rtg.util.StringUtils;
import com.rtg.util.cli.CommandLine;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.io.FileUtils;

/**
 * Output processor for protein matching phase, does multi-core alignment in thread safe manner
 * and produces single output file.
 *
 */
public abstract class ProteinOutputProcessor implements OutputProcessor {
  /** file to store alignments */
  public static final String TABULAR_ALIGNMENTS = "alignments.tsv";

  /** Header line prefix used to identify the line containing the output format version */
  static final String MAPX_OUTPUT_VERSION_HEADER = "#MAPX output";

  /** Header line prefix used to identify the line containing read SDF id */
  static final String MAPX_READ_SDF_ID_HEADER = "#READ-SDF-ID";

  /** Header line prefix used to identify the line containing template SDF id */
  static final String MAPX_TEMPLATE_SDF_ID_HEADER = "#TEMPLATE-SDF-ID";

  /** Current output format version being produced */
  static final String MAPX_OUTPUT_VERSION = "v2.0";

  /** Column header for read-id based alignments */
  static final String HEADER_COL_NAME_READID = "read-id";

  /** Column header for read-name based alignments */
  static final String HEADER_COL_NAME_READNAME = "read-name";

  /** Column headers for read-id based alignments */
  private static final String[] HEADER_COL_NAMES = {
    "template-name",
    "frame",
    HEADER_COL_NAME_READID,
    "template-start",
    "template-end",
    "template-length",
    "read-start",
    "read-end",
    "read-length",
    "template-protein",
    "read-protein",
    "alignment",
    "identical",
    "%identical",
    "positive",
    "%positive",
    "mismatches",
    "raw-score",
    "bit-score",
    "e-score"
  };


  /** specific column indexes */
  private static final int READ_ID_COL = 2;
  private static final int TEMPLATE_PROTEIN_COL = 9;
  private static final int ALIGNMENT_COL = 11;

  /** Column headers for read-id based unmapped file */
  static final String UNMAPPED_HEADER = "#read-id\treason-unmapped";
  /** Column headers for read-name based unmapped file */
  static final String UNMAPPED_HEADER_READ_NAMES = "#read-name\treason-unmapped";

  static final int[] INTERNAL_ENCODED_FRAME_TO_NATURAL_FRAME = {1, 2, 3, -1, -2, -3};
  // FRAMES_MAPPING[frame + 3] gives appropriate TranslatedFrame
  static final Frame[] FRAMES_MAPPING = {
    TranslatedFrame.REVERSE3,
    TranslatedFrame.REVERSE2,
    TranslatedFrame.REVERSE1,
    null,
    TranslatedFrame.FORWARD1,
    TranslatedFrame.FORWARD2,
    TranslatedFrame.FORWARD3};
  private static final Frame[] FRAMES = TranslatedFrame.values();
  //static final String ALIGNMENT_OUT = TABULAR_ALIGNMENTS;
  private static final String ALIGNMENT_SUB_OUTPUT = "alignmentSubOutput_";
  private static final byte[] LS = StringUtils.LS.getBytes();
  static final String UNMAPPED_FILE = "unmapped.tsv";
  private static final int SHIFT_LIMIT = 5;

  // hit pre-filter settings
  private final ProteinBitScorer mPreFilter;
  private final int mPreFilterAlgorithm;
  private final int mPreFilterMinScorePercentage;
  private final int mPreFilterMinOverlapPercentage;

  protected final NgsParams mParams;
  private final SequencesReader mTemplate;
  private final SequencesReader mRead;
  private final ProteinScoringMatrix mProteinScoringMatrix;
  private final BidirectionalEditDistance mProteinEditDistance;
  private final OutputStream mOut;
  private final IntegerOrPercentage mThreshold;
  protected ArrayList<ProteinOutputProcessor> mChildren;
  protected final File mOutFile;
  protected final ProteinOutputProcessor mMaster;
  private final SharedProteinResources mSharedResources;

  private final byte[] mReadWorkspace;
  private long mCurrentTemplateId = -1;
  private byte[] mCurrentTemplate = null;
  protected final SharedStatusCollector mSharedStatusCollector;
  protected final MapStatistics mStatistics;
  private final boolean mEnableReadCache;
  private final byte[] mProteinWorkspace; // = new byte[100]; // Used when not using the read cache

  private long mNumberOfAlignments = 0;
  private long mNumberCacheHits = 0;
  private long mAlignmentsRetained = 0;
  private final long[] mAlignmentShiftCount = new long[21];
  private int mFailedAlignmentThresholdCount = 0;
  private int mFailedPercentIdThresholdCount = 0;
  private int mFailedDueToEscorCount = 0;
  private int mFailedDueToBitScoreCount = 0;
  private int mSkippedDueToFastIdentityFilter = 0;
  private int mSkippedDueToStartLocation = 0;
  private int mAlignmentsDone = 0;
  private int mAlignmentsRepeated = 0;

  private int mLastReadLength = -1;
  private int mLastReadLengthMaxShiftValue = -1;


  /**
   * Create a new master {@link ProteinOutputProcessor}
   * @param params {@link NgsParams} object
   * @param statistics collector of statistics.
   * @throws IOException if error
   */
  protected ProteinOutputProcessor(final NgsParams params, MapStatistics statistics) throws IOException {
    this(params, null, 0, false, null, statistics);
  }

  /**
   * Create a new {@link ProteinOutputProcessor}, internal constructor.
   * @param params {@link NgsParams} object
   * @param master when creating child processors, pass their master in here.  Else null.
   * @param childId the id of this child, starting from 0
   * @param stfu tells non-master output processors to not bother setting up output files.
   * @param collector tells non-master to keep output status.  Null for master processors.
   * @param statistics collector of statistics.
   * @throws IOException if error
   */
  protected ProteinOutputProcessor(final NgsParams params, final ProteinOutputProcessor master, final int childId, final boolean stfu, final SharedStatusCollector collector, MapStatistics statistics) throws IOException {
    mPreFilterMinScorePercentage = params.outputParams().filter().preFilterMinScore();
    mPreFilterMinOverlapPercentage = params.outputParams().filter().preFilterMinOverlap();
    mPreFilterAlgorithm = params.outputParams().filter().preFilterAlgorithm();
    mPreFilter = new ProteinBitScorer(Math.abs(mPreFilterAlgorithm));
    mParams = params;
    mTemplate = mParams.searchParams().reader().copy();
    mRead = mParams.buildFirstParams().reader().copy();
    //assert mTemplate instanceof CompressedMemorySequencesReader;
    //assert mRead instanceof CompressedMemorySequencesReader;
    mReadWorkspace = new byte[(int) mRead.maxLength()];
    mEnableReadCache = mParams.enableProteinReadCache();
    mProteinScoringMatrix = mParams.proteinScoringMatrix();
    mProteinEditDistance = EditDistanceFactory.createProteinEditDistance(mProteinScoringMatrix);
    mSharedResources = new SharedProteinResources(mProteinScoringMatrix, mTemplate, mRead, mParams.outputParams().outputReadNames());
    mThreshold = mParams.outputParams().filter().matedMaxMismatches();
    mStatistics = statistics;
    if (master == null) {
      // this is the master
      mMaster = this;
      if (collector != null) {
        throw new IllegalArgumentException("collector can not be preset");
      }
      final boolean compressed = params.outputParams().isCompressOutput();
      final File outDir = mParams.outputParams().directory();
      mOutFile = new File(outDir, compressed ? TABULAR_ALIGNMENTS + FileUtils.GZ_SUFFIX : TABULAR_ALIGNMENTS);
      mOut = FileUtils.createOutputStream(mOutFile);
      writeAlignmentHeader(mOut);
      mChildren = new ArrayList<>();
      assert mRead.numberSequences() <= Integer.MAX_VALUE;
      mSharedStatusCollector = new SharedStatusCollector((int) mRead.numberSequences(), mStatistics);
    } else {
      // this is a child, so remember who its master is.
      mMaster = master;
      mChildren = null;
      if (collector == null) {
        throw new IllegalArgumentException("children needs shared read status collector");
      }
      mSharedStatusCollector = collector;
      if (!stfu) {
        final String name = ALIGNMENT_SUB_OUTPUT + childId;
        final File tempDir = mParams.outputParams().tempFilesDirectory();
        if (!tempDir.exists() && !tempDir.mkdirs()) {
          throw new IOException("Could not create temporary directory: " + tempDir.getPath());
        }
        mOutFile = params.outputParams().resultStreamHandler().tempFile(name + (mParams.outputParams().isCompressOutput() ? FileUtils.GZ_SUFFIX : ""));
        mOut = FileUtils.createOutputStream(mOutFile);
        if (childId == 0) {
          writeAlignmentHeader(mOut);
        }
      } else {
        mOutFile = null;
        mOut = null;
      }
    }
    mProteinWorkspace = new byte[(int) params.buildFirstParams().reader().maxLength()];
  }

  protected NgsParams getParams() {
    return mParams;
  }

  /**
   * Return true if the current result could be retained.  This allows for immediate
   * filtering of alignments prior to <code>ProteinAlignmentResult</code> computation.
   * @param res packed protein result
   * @param readId read id
   * @param templateId template id
   * @param readAndFrame read id and frame
   * @param plen length in protein space
   * @return false if result should be discarded
   * @throws java.io.IOException if an I/O error occurs
   */
    protected boolean retainResult(final int[] res, final int readId, final int templateId, final int readAndFrame, final int plen) throws IOException {
       final int alignmentScore = ActionsHelper.alignmentScore(res);
//    boolean retain = true;
    if (alignmentScore > mThreshold.getValue(plen)) {
      mSharedStatusCollector.setStatus(readId, SharedStatusCollector.EXCEEDS_ALIGNMENT_THRESHOLD);
      ++mFailedAlignmentThresholdCount;
      return false;
    }
    if (ScoringHelper.percentId(res) < mParams.outputParams().filter().minIdentity()) {
      mSharedStatusCollector.setStatus(readId, SharedStatusCollector.EXCEEDS_PERCENT_ID_THRESHOLD);
      ++mFailedPercentIdThresholdCount;
      return false;
    }
    if (ScoringHelper.computeEScore(alignmentScore, ActionsHelper.actionsCount(res), mSharedResources.totalTemplateLength(), mProteinScoringMatrix) > mParams.outputParams().filter().maxEScore()) {
      mSharedStatusCollector.setStatus(readId, SharedStatusCollector.EXCEEDS_E_SCORE_THRESHOLD);
      ++mFailedDueToEscorCount;
      return false;
    }
    if (ScoringHelper.computeBitScore(alignmentScore, mProteinScoringMatrix) < mParams.outputParams().filter().minBitScore()) {
      mSharedStatusCollector.setStatus(readId, SharedStatusCollector.BELOW_BIT_SCORE_THRESHOLD);
      ++mFailedDueToBitScoreCount;
      return false;
    }
   return true;
  }

  @Override
  public void process(final long templateId, final String fr, final int rr, final int chunkStart, final int score, final int scoreIndel) throws IOException {
    if (mCurrentTemplateId != templateId) {
      nextTemplateId(templateId);
    }
    final int chunkId = rr / ((int) mRead.numberSequences() * FRAMES.length);
    final int r = rr - chunkId * (int) mRead.numberSequences() * FRAMES.length;
    final int readId = r / FRAMES.length;

    final int genomeFrame = INTERNAL_ENCODED_FRAME_TO_NATURAL_FRAME[r % FRAMES.length];
    final Frame frames = FRAMES_MAPPING[genomeFrame + 3];
    //Diagnostic.developerLog("readId: " + readId);
    final int plen;
    byte[] readProtein = mEnableReadCache ? mSharedStatusCollector.getReadProtein(r) : null;
    if (readProtein == null) {
      // get read, convert to protein, put it in the cache
      final int rlen = mRead.read(readId, mReadWorkspace);
      plen = (rlen - Math.abs(genomeFrame) + 1) / 3;
      readProtein = mEnableReadCache ? new byte[plen] : mProteinWorkspace;
      for (int j = 0, i = 0; j < plen; ++j, i += 3) {
        readProtein[j] = frames.code(mReadWorkspace, rlen, i);
      }
      if (mEnableReadCache) {
        mSharedStatusCollector.putReadProtein(r, readProtein);
      }
    } else {
      plen = readProtein.length;
      ++mNumberCacheHits;
    }


    final int chunkToPosition = ProteinReadIndexer.chunkToPosition(chunkId, mRead.length(readId) / 3, mParams.mapXMetaChunkSize(), mParams.mapXMetaChunkOverlap());
    final int chunkReadStart;
    final int start;
    if (frames.isForward()) {
      start = chunkStart - chunkToPosition;
      chunkReadStart = chunkToPosition;
    } else {
      int chunkEnd = chunkToPosition + mParams.mapXMetaChunkSize();
      if (chunkEnd > plen) {
        chunkEnd = plen;
      }
      start = chunkStart - (plen - chunkEnd);
      chunkReadStart = Math.max(plen - chunkToPosition - mParams.mapXMetaChunkSize(), 0);
    }

    ++mNumberOfAlignments;
    final int minOverlap = mPreFilterMinOverlapPercentage * plen / 100;
    if (start + plen < minOverlap || mCurrentTemplate.length - start < minOverlap) {
      ++mSkippedDueToStartLocation;
      return;
    }
    final int fastScoreReadLength = Math.min(plen, mParams.mapXMetaChunkSize());
    final int minScore = mPreFilterMinScorePercentage * fastScoreReadLength / 100;
    final int[] bitScores = mPreFilter.calculateFastScore(readProtein, chunkReadStart, fastScoreReadLength, mCurrentTemplate, chunkStart);
    final boolean skipIt;
    if (mPreFilterAlgorithm < 0) {
      // We look for the highest peak, like 2 1 14 3 2
      // We don't care about its exact position, since the aligner will sort that out.
      final int width = bitScores.length - 1;
      int sum = 0;
      for (int i = 0; i < width; ++i) {
        sum += bitScores[i];
      }
      // now search for the highest peak
      int maxPeak = 0;
      for (int i = 0; i < width; ++i) {
        // peak = sum{j:0..width-1, excluding i} of bitScores[i] - bitScores[j]
        final int peak = bitScores[i] * width - sum;
        if (peak > maxPeak) {
          maxPeak = peak;
        }
      }
      // we divide maxPeak by (width / 2) so that minScore is independent of width.
      skipIt = maxPeak * 2 / width < minScore;
      //Diagnostic.developerLog("maxPeak: " + maxPeak + " width: " + width + " calc: " + (maxPeak * 2 / width) + " minScore: " + minScore);
    } else {
      // the old bitOR algorithm
      skipIt = bitScores[bitScores.length - 1] < minScore;
      //Diagnostic.developerLog("bitScores: " + bitScores[bitScores.length - 1] + " minScore: " + minScore);
    }
    if (skipIt) {
      ++mSkippedDueToFastIdentityFilter;
      return;
    }
    if (mLastReadLength != plen) {
      mLastReadLength = plen;
      mLastReadLengthMaxShiftValue = calculateProteinMaxShift(plen);
    }
    ++mAlignmentsDone;
    int[] res = mProteinEditDistance.calculateEditDistance(readProtein, plen, mCurrentTemplate, start, false, Integer.MAX_VALUE, mLastReadLengthMaxShiftValue, true);

    int newZeroBasedTemplateStart = ActionsHelper.zeroBasedTemplateStart(res);
    if (Math.abs(newZeroBasedTemplateStart - start) > SHIFT_LIMIT) {
      ++mAlignmentsRepeated;
      res = mProteinEditDistance.calculateEditDistance(readProtein, plen, mCurrentTemplate, newZeroBasedTemplateStart, false, Integer.MAX_VALUE, mLastReadLengthMaxShiftValue, true);
      newZeroBasedTemplateStart = ActionsHelper.zeroBasedTemplateStart(res);
    }
    final int shift = newZeroBasedTemplateStart - start;

    if (retainResult(res, readId, (int) templateId, r, plen)) {
      ++mAlignmentsRetained;

      final long idOffset = Math.max(mParams.buildFirstParams().readerRestriction().getStart(), 0);
      writeResult(new ProteinAlignmentResult(mSharedResources, (int) templateId, r, res, idOffset, mParams.outputParams().outputProteinSequences()));
      // and record a histogram of the shift distances
      if (Math.abs(shift) > mAlignmentShiftCount.length / 2) {
        Diagnostic.developerLog("alignment for read " + r + ", template " + templateId + ":" + start + " shifted by " + shift);
      } else {
        mAlignmentShiftCount[shift + mAlignmentShiftCount.length / 2]++;
      }
    }
  }

  /**
   * This is done based on information available at this time, for 100 long DNA reads (32/33 long AA reads)
   * biggest insertion found is 18AA approx. 56% of read length, the current value will be calculated using
   * <code>readlen * 0.6</code>
   * @param plen length of the read in protein space
   * @return allowed max shift
   */
  private int calculateProteinMaxShift(int plen) {
    return (int) (plen * 0.6 + 0.5);
  }

  /**
   * Switch to another template
   * @param templateId the id of the next template sequence
   * @throws java.io.IOException if an I/O error occurs
   */
  protected void nextTemplateId(long templateId) throws IOException {
    mCurrentTemplateId = templateId;
    mCurrentTemplate = mTemplate.read(templateId);
  }

  void writeResult(final ProteinAlignmentResult res) throws IOException {
    mSharedStatusCollector.setStatus(res.readId(), SharedStatusCollector.RESULT_WRITTEN);
    res.write(mOut);
    mOut.write(LS);
  }

  @Override
  public void close() throws IOException {
    if (mOut != null) {
      mOut.close();
    }
    Diagnostic.developerLog("Fast identity percentage : " + mPreFilterMinScorePercentage + "%");
    Diagnostic.developerLog("Overhang identity percentage : " + mPreFilterMinOverlapPercentage + "%");
    Diagnostic.developerLog("Total alignments : " + mNumberOfAlignments);
    Diagnostic.developerLog("Read cache hits  : " + mNumberCacheHits);
    Diagnostic.developerLog("Alignments skipped due to offset start   : " + mSkippedDueToStartLocation);
    Diagnostic.developerLog("Alignments skipped due to fast identity  : " + mSkippedDueToFastIdentityFilter);
    Diagnostic.developerLog("Alignments done      : " + mAlignmentsDone);
    Diagnostic.developerLog("Alignments done twice: " + mAlignmentsRepeated);
    Diagnostic.developerLog("Alignments failed due to alignment score : " + mFailedAlignmentThresholdCount);
    Diagnostic.developerLog("Alignments failed due to percent ID      : " + mFailedPercentIdThresholdCount);
    Diagnostic.developerLog("Alignments failed due to E score         : " + mFailedDueToEscorCount);
    Diagnostic.developerLog("Alignments failed due to bit score       : " + mFailedDueToBitScoreCount);
    Diagnostic.developerLog("Alignments retained : " + mAlignmentsRetained);

    for (int i = 0; i < mAlignmentShiftCount.length; ++i) {
      final int shift = i - mAlignmentShiftCount.length / 2;
      Diagnostic.developerLog("Alignments shifted by " + shift + " : " + mAlignmentShiftCount[i]);
    }
    mProteinEditDistance.logStats();
  }


  private void writeAlignmentHeader(OutputStream out) throws IOException {
    writeCommonHeader(out);
    final StringBuilder sb = new StringBuilder();
    sb.append('#').append(HEADER_COL_NAMES[0]);
    for (int i = 1; i < HEADER_COL_NAMES.length; ++i) {
      if (i == READ_ID_COL && mParams.outputParams().outputReadNames()) {
        sb.append('\t').append(HEADER_COL_NAME_READNAME);
      } else if (i < TEMPLATE_PROTEIN_COL
          || i > ALIGNMENT_COL
          || mParams.outputParams().outputProteinSequences()) {
        sb.append('\t').append(HEADER_COL_NAMES[i]);
      }
    }
    out.write(sb.toString().getBytes());
    out.write(LS);
  }

  private void writeUnmappedHeader(OutputStream out) throws IOException {
    writeCommonHeader(out);
    out.write((mParams.outputParams().outputReadNames() ? UNMAPPED_HEADER_READ_NAMES : UNMAPPED_HEADER).getBytes());
    out.write(LS);
  }

  private void writeCommonHeader(OutputStream out) throws IOException {
    out.write(("#Version\t" + Environment.getVersion()).getBytes());
    out.write(LS);
    out.write((MAPX_OUTPUT_VERSION_HEADER + "\t" + MAPX_OUTPUT_VERSION).getBytes());
    out.write(LS);
    if (CommandLine.getCommandLine() != null) {
      out.write(("#CL\t" + CommandLine.getCommandLine()).getBytes());
      out.write(LS);
    }
    out.write(("#RUN-ID\t" + CommandLine.getRunId() + StringUtils.LS).getBytes());
    if (mTemplate.getSdfId().available()) {
      out.write((MAPX_TEMPLATE_SDF_ID_HEADER + "\t" + mTemplate.getSdfId()).getBytes());
      out.write(LS);
    }
    if (mRead.getSdfId().available()) {
      out.write((MAPX_READ_SDF_ID_HEADER + "\t" + mRead.getSdfId()).getBytes());
      out.write(LS);
    }
  }

  void writeUnmapped() throws IOException {
    if (mParams.outputParams().outputUnmapped()) {
      final boolean compressed = mParams.outputParams().isCompressOutput();
      final File outDir = mParams.outputParams().directory();
      try (OutputStream unmapped = FileUtils.createOutputStream(new File(outDir, compressed ? UNMAPPED_FILE + FileUtils.GZ_SUFFIX : UNMAPPED_FILE))) {
        writeUnmappedHeader(unmapped);
        final long idOffset = Math.max(mParams.buildFirstParams().readerRestriction().getStart(), 0);
        mSharedStatusCollector.writeUnmapped(unmapped, mParams.outputParams().outputReadNames() ? mParams.buildFirstParams().reader().names() : null, idOffset);
      }
    }
  }

  protected void closeChildren() throws IOException {
    for (final ProteinOutputProcessor child : mChildren) {
      child.finish();
      child.close();
    }
  }

}
