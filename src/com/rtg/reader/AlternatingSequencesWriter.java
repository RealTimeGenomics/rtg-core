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
import java.util.Collection;

import com.rtg.util.intervals.LongRange;
import com.rtg.util.cli.CommandLine;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.NoTalkbackSlimException;
import com.rtg.util.diagnostic.WarningType;


/**
 * Writes sequences into two SDFs, alternating between left and right.
 * Used in CG data formatting.
 */
public class AlternatingSequencesWriter extends SequencesWriter {

  /**
   * Creates an alternating writer for processing sequences from provided data source.
   * @param source Source of sequences
   * @param outputDir Destination of output files
   * @param sizeLimit Maximum size for output files
   * @param type preread type
   * @param compressed whether <code>SDF</code> should be compressed
   * @param trimQualityThreshold quality threshold to trim reads after, or null for no trimming
   */
  public AlternatingSequencesWriter(SequenceDataSource source, File outputDir, long sizeLimit, PrereadType type, boolean compressed, Integer trimQualityThreshold) {
    super(source, outputDir, sizeLimit, null, type, compressed, trimQualityThreshold);
  }
  /**
   * Creates an alternating writer for processing sequences from provided data source.
   * @param source Source of sequences
   * @param outputDir Destination of output files
   * @param sizeLimit Maximum size for output files
   * @param namesToExclude Names to be excluded from written result (may be null)
   * @param type preread type
   * @param compressed whether <code>SDF</code> should be compressed
   * @param trimQualityThreshold quality threshold to trim reads after, or null for no trimming
   */
  public AlternatingSequencesWriter(SequenceDataSource source, File outputDir, long sizeLimit, final Collection<String> namesToExclude, PrereadType type, boolean compressed, Integer trimQualityThreshold) {
    super(source, outputDir, sizeLimit, namesToExclude, type, compressed, trimQualityThreshold);
  }

  /**
   * Constructor for use with in memory sequence sinks
   * @param source Source of sequences
   * @param namesToExclude Names to be excluded from written result (may be null)
   * @param type preread type
   * @param compressed whether <code>SDF</code> should be compressed
   */
  public AlternatingSequencesWriter(final SequenceDataSource source, final Collection<String> namesToExclude, final PrereadType type, boolean compressed) {
    super(source, namesToExclude, type, compressed);
  }

  private void processSequences(AbstractSdfWriter sdfWriterLeft, AbstractSdfWriter sdfWriterRight, LongRange region) throws IOException {
    try {
      final boolean checkLimit = region != null && region != LongRange.NONE && region.getEnd() != LongRange.MISSING;
      final long seqLimit;
      if (checkLimit) {
        seqLimit = region.getEnd();
      } else {
        seqLimit = -1;
      }
      sdfWriterLeft.setPrereadArm(PrereadArm.LEFT);
      sdfWriterLeft.setSdfId(mSdfId);
      sdfWriterLeft.setComment(mComment);
      sdfWriterLeft.setReadGroup(mReadGroup);
      sdfWriterLeft.setCommandLine(CommandLine.getCommandLine());
      sdfWriterRight.setPrereadArm(PrereadArm.RIGHT);
      sdfWriterRight.setSdfId(mSdfId);
      sdfWriterRight.setComment(mComment);
      sdfWriterRight.setReadGroup(mReadGroup);
      sdfWriterRight.setCommandLine(CommandLine.getCommandLine());

      String name = "";
      while (mDataSource.nextSequence()) {
        //this is to short circuit processing
        if (checkLimit && sdfWriterLeft.getNumberOfSequences() >= seqLimit) {
          break;
        }
        if (mCheckDuplicateNames) {
          if (name.equals(mDataSource.name())) {
            throw new NoTalkbackSlimException("More than one read-pair with the same name in input.");
          }
          name = mDataSource.name();
        }
        processSingleSequence(sdfWriterLeft);
        if (!mDataSource.nextSequence()) {
          throw new NoTalkbackSlimException("Input source had uneven number of records.");
        }
        processSingleSequence(sdfWriterRight);
      }

      if (mDataSource.getWarningCount() > FastaSequenceDataSource.NUMBER_OF_TIDE_WARNINGS) {
        Diagnostic.warning(""); //Ugly way to separate the warning.
        Diagnostic.warning(WarningType.NUMBER_OF_BAD_TIDE, Long.toString(mDataSource.getWarningCount()));
      }
    } finally {
      try {
        sdfWriterLeft.close();
      } finally {
        sdfWriterRight.close();
      }
      calculateStats(sdfWriterLeft, sdfWriterRight);
    }
  }

  /**
   * Processes alternating paired input sequences into two in memory readers
   * @param sourceFile file that was the source of this data
   * @param includeQuality true if including quality data in output, false otherwise
   * @param names storage for read names
   * @param suffixes storage for read name suffixes
   * @param region restriction returned reader to given range of reads
   * @return the reader from reading the sequence data
   * @throws IOException If an I/O error occurs
   */
  public CompressedMemorySequencesReader[] processSequencesInMemoryPaired(File sourceFile, boolean includeQuality, SimplePrereadNames names, SimplePrereadNames suffixes, LongRange region) throws IOException {
    try {
      final CompressedMemorySequencesWriter sdfWriterLeft = new CompressedMemorySequencesWriter(sourceFile, mPrereadType, mDataSource.hasQualityData() && includeQuality, names, suffixes, true, mDataSource.type(), region);
      final RightSimplePrereadNames rNames = names == null ? null : new RightSimplePrereadNames(names);
      final RightSimplePrereadNames rSuffixes = names == null ? null : new RightSimplePrereadNames(suffixes);
      final CompressedMemorySequencesWriter sdfWriterRight = new CompressedMemorySequencesWriter(sourceFile, mPrereadType, mDataSource.hasQualityData() && includeQuality, rNames, rSuffixes, true, mDataSource.type(), region);
      processSequences(sdfWriterLeft, sdfWriterRight, region);
      return new CompressedMemorySequencesReader[] {sdfWriterLeft.getReader(), sdfWriterRight.getReader()};
    } finally {
      mDataSource.close();
    }
  }


  @Override
  public void processSequences(boolean includeQuality, boolean includeNames) throws IOException {
    final File outputDirLeft = new File(mOutputDir, "left");
    final File outputDirRight = new File(mOutputDir, "right");
    try {
      final SdfWriter sdfWriterLeft = new SdfWriter(outputDirLeft, mSizeLimit, mPrereadType, mDataSource.hasQualityData() && includeQuality, includeNames, mCompressed, mDataSource.type());
      final SdfWriter sdfWriterRight = new SdfWriter(outputDirRight, mSizeLimit, mPrereadType, mDataSource.hasQualityData() && includeQuality, includeNames, mCompressed, mDataSource.type());
      processSequences(sdfWriterLeft, sdfWriterRight, null);
    } finally {
      mDataSource.close();
    }
  }

  private void calculateStats(AbstractSdfWriter left, AbstractSdfWriter right) {
    mTotalLength += left.getTotalLength() + right.getTotalLength();
    mNumberOfSequences += left.getNumberOfSequences() + right.getNumberOfSequences();
    if (left.getMaxLength() > mMaxLength || right.getMaxLength() > mMaxLength) {
      mMaxLength = Math.max(left.getMaxLength(), right.getMaxLength());
    }
    if (left.getMinLength() < mMinLength || right.getMinLength() < mMinLength) {
      mMinLength = Math.min(left.getMinLength(), right.getMinLength());
    }
  }

}
