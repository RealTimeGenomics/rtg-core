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
package com.rtg.reader;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;

import com.rtg.mode.DNAFastaSymbolTable;
import com.rtg.util.Constants;
import com.rtg.util.cli.CommandLine;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.NoTalkbackSlimException;
import com.rtg.util.diagnostic.WarningType;
import com.rtg.util.intervals.LongRange;
import com.rtg.util.io.FileUtils;

/**
 * Takes a sequence data source, and writes sequences out in the preread data format.
 * Assumes the destination directory is essentially empty and we can do whatever we want in it.
 */
public class SequencesWriter {

  final SequenceDataSource mDataSource;
  final File mOutputDir;
  long mSizeLimit = SdfWriter.DEFAULT_SIZE_LIMIT;

  //Stats
//  private static final int MAX_WARNINGS = 5;
//  private final boolean mSkipEmptySequences = false; //Boolean.valueOf(System.getProperty("skip-empty-sequences", "false")); // If true, don't pass empty sequences to sdf writer
  private int mEmptySequenceCount = 0;

  private long mExcludedSequencesCount = 0;
  private long mExcludedResiduesCount;
  final PrereadType mPrereadType;
  private PrereadArm mPrereadArm = PrereadArm.UNKNOWN;

  SdfId mSdfId;

  String mComment = null;
  String mReadGroup;

  private final Collection<String> mNamesToExclude;

  /** Set to true when input source is sorted on read name and you want to check no consecutive reads or pairs of reads have the same name */
  protected boolean mCheckDuplicateNames;

  long mTotalLength = 0;
  long mMinLength = Long.MAX_VALUE;
  long mMaxLength = Long.MIN_VALUE;
  long mNumberOfSequences = 0;
  final boolean mCompressed;

  final ReadTrimmer mReadTrimmer;

  /**
   * Creates a writer for processing sequences from provided data source.
   * @param source Source of sequences
   * @param outputDir Destination of output files
   * @param sizeLimit Maximum size for output files
   * @param namesToExclude Names to be excluded from written result (may be null)
   * @param type preread type
   * @param trimQualityThreshold quality threshold to trim reads after, or null for no trimming
   * @param compressed whether <code>SDF</code> should be compressed
   */
  public SequencesWriter(SequenceDataSource source, File outputDir, long sizeLimit, Collection<String> namesToExclude, PrereadType type, boolean compressed, Integer trimQualityThreshold) {
    mDataSource = source;
    mOutputDir = outputDir;
    mPrereadType = type;
    mSdfId = new SdfId();
    if (mDataSource == null) {
      throw new NullPointerException("Provided data source was null.");
    }
    FileUtils.ensureOutputDirectory(mOutputDir);
    mSizeLimit = sizeLimit;
    mNamesToExclude = namesToExclude == null ? new ArrayList<String>() : namesToExclude;
    mCompressed = compressed;
    if (trimQualityThreshold != null) {
      if (System.getProperty("gillette-window-size") != null) {
        final int windowSize = Integer.parseInt(System.getProperty("gillette-window-size"));
        mReadTrimmer = new GilletteReadTrimmer(windowSize, trimQualityThreshold);
        Diagnostic.userLog("Performing trimming with Gillette, threshold " + trimQualityThreshold + " window size " + windowSize);
      } else {
        mReadTrimmer = new BestSumReadTrimmer(trimQualityThreshold);
        Diagnostic.userLog("Performing trimming with BestSum, threshold " + trimQualityThreshold);
      }
    } else {
      mReadTrimmer = null;
    }
  }

  /**
   * Constructor for use with in memory sequence sinks
   * @param source Source of sequences
   * @param namesToExclude Names to be excluded from written result (may be null)
   * @param type preread type
   * @param compressed whether <code>SDF</code> should be compressed
   */
  public SequencesWriter(SequenceDataSource source, Collection<String> namesToExclude, PrereadType type, boolean compressed) {
    mDataSource = source;
    mPrereadType = type;
    mSdfId = new SdfId();
    if (mDataSource == null) {
      throw new NullPointerException("Provided data source was null.");
    }
    mOutputDir = null;
    mSizeLimit = 0;
    mNamesToExclude = namesToExclude == null ? new ArrayList<String>() : namesToExclude;
    mCompressed = compressed;
    mReadTrimmer = null;
  }

  /**
   * Creates a writer for processing sequences from provided data source.
   * @param source Source of sequences
   * @param outputDir Destination of output files
   * @param sizeLimit Maximum size for output files
   * @param type preread type
   * @param compressed whether <code>SDF</code> should be compressed
   */
  public SequencesWriter(SequenceDataSource source, File outputDir, long sizeLimit, PrereadType type, boolean compressed) {
    this(source, outputDir, sizeLimit, null, type, compressed, null);
  }

  /**
   * Set to true when input source is sorted on read name and you want to check no consecutive reads or pairs of reads have the same name.
   * @param val true when input source is sorted on read name and you want to check no consecutive reads or pairs of reads have the same name, false otherwise.
   */
  public void setCheckDuplicateNames(boolean val) {
    mCheckDuplicateNames = val;
  }

  /**
   * Returns the size limit
   * @return the size limit of the file
   */
  public long getSizeLimit() {
    return mSizeLimit;
  }

//  private void noSequenceWarning(String label) {
//    if (mEmptySequenceCount < MAX_WARNINGS) {
//      Diagnostic.warning(WarningType.NO_SEQUENCE, label);
//    }
//    mEmptySequenceCount++;
//    if (mEmptySequenceCount == MAX_WARNINGS) {
//      Diagnostic.warning("Subsequent warnings of this type will not be shown.");
//    }
//  }

  /**
   * Writes a single sequence.
   * @throws IOException if an I/O error occurs.
   */
  void processSingleSequence(final AbstractSdfWriter mSdfWriter) throws IOException {
    final String label = mDataSource.name();

    //filtering, not related to format of data
    for (final String s : mNamesToExclude) {
      if (label.contains(s)) {
        mExcludedSequencesCount++;
        mExcludedResiduesCount += mDataSource.currentLength();
        // Found a sequence that is supposed to be suppressed.
        //while (mDataSource.read() != null) {
        //  ;
        //}
        return;
      }
    }

    final int length = mReadTrimmer != null && mDataSource.hasQualityData() ? mReadTrimmer.getTrimPosition(mDataSource.qualityData(), mDataSource.currentLength()) : mDataSource.currentLength();
    mExcludedResiduesCount += mDataSource.currentLength() - length;
//    if (length == 0 && mSkipEmptySequences) {
//      noSequenceWarning(label);
//    } else {
      mSdfWriter.startSequence(label);
      mSdfWriter.write(mDataSource.sequenceData(), mDataSource.qualityData(), length);  //may write 0 length
      mSdfWriter.endSequence();
//    }
  }

  /**
   * Processes input sequences into output files
   * @throws IOException If an I/O error occurs
   */
  public void processSequences() throws IOException {
    processSequences(true, true);
  }

  private void processSequences(AbstractSdfWriter inSdfWriter, LongRange region) throws IOException {
    try (AbstractSdfWriter sdfWriter = inSdfWriter;
         SequenceDataSource dataSource = mDataSource) {
      final boolean checkLimit = region != null && region != LongRange.NONE && region.getEnd() != LongRange.MISSING;
      final long seqLimit;
      if (checkLimit) {
        seqLimit = region.getEnd();
      } else {
        seqLimit = -1;
      }
      sdfWriter.setPrereadArm(mPrereadArm);
      sdfWriter.setSdfId(mSdfId);
      sdfWriter.setComment(mComment);
      sdfWriter.setReadGroup(mReadGroup);
      sdfWriter.setCommandLine(CommandLine.getCommandLine());
      String name = "";
      while (dataSource.nextSequence()) {
        //this is to short circuit processing
        if (checkLimit && sdfWriter.getNumberOfSequences() >= seqLimit) {
          break;
        }
        if (mCheckDuplicateNames) {
          if (name.equals(dataSource.name())) {
            throw new NoTalkbackSlimException("More than one read with the same name in input.");
          }
          name = dataSource.name();
        }
        processSingleSequence(sdfWriter);
      }
      if (dataSource.getWarningCount() > FastaSequenceDataSource.NUMBER_OF_TIDE_WARNINGS) {
        Diagnostic.warning(""); //Ugly way to separate the warning.
        Diagnostic.warning(WarningType.NUMBER_OF_BAD_TIDE, Long.toString(dataSource.getWarningCount()));
      }
      if (mEmptySequenceCount > 0) {
        Diagnostic.warning(""); //Ugly way to separate the warning.
        Diagnostic.warning("There were " + mEmptySequenceCount + " sequences with no data.");
      }
    } finally {
      mTotalLength += inSdfWriter.getTotalLength();
      mNumberOfSequences += inSdfWriter.getNumberOfSequences();
      mMaxLength = inSdfWriter.getMaxLength();
      mMinLength = inSdfWriter.getMinLength();
    }
  }

  /**
   * Processes input sequences into memory
   * @param sourceFile if the data originates in a file specify it here otherwise null.
   * @param includeQuality true if including quality data in output, false otherwise
   * @param names read name storage
   * @param suffixes read name suffix storage
   * @param region restriction returned reader to given range of reads
   * @return the reader from reading the sequence data
   * @throws IOException If an I/O error occurs
   */
  public CompressedMemorySequencesReader processSequencesInMemory(File sourceFile, boolean includeQuality, SimplePrereadNames names, SimplePrereadNames suffixes, LongRange region) throws IOException {
    final CompressedMemorySequencesWriter sdfWriter = new CompressedMemorySequencesWriter(sourceFile, mPrereadType, mDataSource.hasQualityData() && includeQuality, names, suffixes, true, mDataSource.type(), region);
    processSequences(sdfWriter, region);
    return sdfWriter.getReader();
  }

  /**
   * Processes input sequences into output files
   * @param includeQuality true if including quality data in output, false otherwise
   * @param includeNames true if including name data in output, false otherwise
   * @throws IOException If an I/O error occurs
   */
  public void processSequences(boolean includeQuality, boolean includeNames) throws IOException {
    final SdfWriter sdfWriter = new SdfWriter(mOutputDir, mSizeLimit, mPrereadType, mDataSource.hasQualityData() && includeQuality, includeNames, mCompressed, mDataSource.type());
    processSequences(sdfWriter, null);
  }

  /**
   * @param prereadArm the preread arm to set
   */
  public final void setPrereadArm(final PrereadArm prereadArm) {
    mPrereadArm = prereadArm;
  }

  /**
   * @param id set the SDF-ID for this sequence
   */
  public void setSdfId(final SdfId id) {
    mSdfId = id;
  }

  /**
   * @param comment the comment for the SDF file
   */
  public void setComment(final String comment) {
    mComment = comment;
  }

  /**
   * @param readGroup the sam read group to store in the SDF
   */
  public void setReadGroup(final String readGroup) {
    mReadGroup = readGroup;
  }

  /**
   * @return the total Length
   */
  public final long getTotalLength() {
    return mTotalLength;
  }

  /**
   * @return the maximum length
   */
  public final long getMaxLength() {
    return mMaxLength;
  }

  /**
   * @return the minimum length
   */
  public final long getMinLength() {
    return mMinLength;
  }

  /**
   * @return the number Of Sequences
   */
  public final long getNumberOfSequences() {
    return mNumberOfSequences;
  }

  /**
   * @return the number of Sequences that were excluded
   */
  public final long getNumberOfExcludedSequences() {
    return mExcludedSequencesCount + mEmptySequenceCount;
  }

  /**
   * @return the count of residues from sequences that were excluded
   */
  public final long getExcludedResidueCount() {
    return mExcludedResiduesCount;
  }

  /**
   * @return the SDF Id
   */
  public SdfId getSdfId() {
    return mSdfId;
  }

  /**
   * @param args command line arguments
   * @throws Exception exception
   */
  public static void main(String[] args) throws Exception {
    final long start = System.currentTimeMillis();
    final SequenceDataSource leftds = new FastaSequenceDataSource(Arrays.asList(new File(args[0])), new DNAFastaSymbolTable(), true, null);
    final SequencesWriter writer = new SequencesWriter(leftds, new File(args[1]), Constants.MAX_FILE_SIZE, new ArrayList<String>(), PrereadType.UNKNOWN, true, null);
    // perform the actual work
    writer.processSequences();
    System.err.println("time " + ((System.currentTimeMillis() - start) / 1000) + "s");
    /*
    final long start2 = System.currentTimeMillis();
    final BufferedInputStream sr = FileUtils.createInputStream(new File(args[0]), true);
    final byte[] buf = new byte[1024000];
    while (sr.read(buf) != -1) {
      ;
    }
    System.err.println("sr time " + ((System.currentTimeMillis() - start2) / 1000) + "s");
    */
    System.err.println("done");
  }
}
