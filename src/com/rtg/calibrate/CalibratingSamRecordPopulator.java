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
package com.rtg.calibrate;

import java.io.IOException;
import java.util.Map;

import com.rtg.reader.ReaderUtils;
import com.rtg.reader.SequencesReader;
import com.rtg.sam.SamRecordPopulator;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.NoTalkbackSlimException;
import com.rtg.util.diagnostic.SlimException;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.RuntimeIOException;

/**
 * Populator that passes through SAM records and performs calibration statistic accumulation.
 */
public class CalibratingSamRecordPopulator extends SamRecordPopulator {

  private static final int BUFFER_SIZE = 10000000; // How many bases of template sequence to read at a time

  private final Calibrator mCalibrator;
  private final SequencesReader mReference;
  private final Map<String, Long> mSequenceNameMap;
  private final byte[] mReferenceBytes = new byte[BUFFER_SIZE];

  private String mCurrentName = null;
  private long mCurrentId = 0;
  private int mCurrentLength = 0;

  private int mReferenceStart = 0;
  private int mReferenceEnd = 0;

  /**
   * Populator for SAM records using an initializer. Can be given null to
   * skip initialization.
   * @param calibrator calibrator to use
   * @param reference template reader
   * @param forceInit if true use the default SAM record forced field initializer
   */
  public CalibratingSamRecordPopulator(final Calibrator calibrator, final SequencesReader reference, boolean forceInit) {
    super(forceInit ? SamRecordPopulator.DEFAULT_INITIALIZER : null);
    mCalibrator = calibrator;
    mReference = reference;
    try {
      mSequenceNameMap = ReaderUtils.getSequenceNameMap(mReference);
    } catch (final IOException e) {
      throw new SlimException("Problem getting sequence names from reference");
    }
  }

  @Override
  public SAMRecord populate(final SAMRecord rec) {
    if (!rec.getReadUnmappedFlag()) {
      final String name = rec.getReferenceName();
      try {
        if (!name.equals(mCurrentName)) {
          final Long seqId = mSequenceNameMap.get(name);
          if (seqId == null) {
            //user most likely is calibrating against a different reference to what the mappings came from
            throw new NoTalkbackSlimException("SAM record is aligned to sequence " + name + ", but this is not contained in the supplied reference");
          }
          mCurrentName = name;
          mCurrentId = seqId;
          mCurrentLength = mReference.length(mCurrentId);
          mReferenceStart = 0;
          mReferenceEnd = 0;
        }
        ensureTemplate(rec.getAlignmentStart() - 1, rec.getAlignmentEnd());
      } catch (IOException e) {
        throw new RuntimeIOException(e);
      }
      mCalibrator.processRead(rec);
    }
    return super.populate(rec);
  }

  private void ensureTemplate(int alignmentStart, int alignmentEnd) throws IOException {
    if (alignmentStart < mReferenceStart) {
      throw new IOException("Alignments are not sorted");
    }
    if (mReferenceEnd <= alignmentEnd) {
      mReferenceStart = alignmentStart;
      mReferenceEnd = Math.min(mCurrentLength, Math.max(alignmentEnd, mReferenceStart + mReferenceBytes.length));
      final int length = mReference.read(mCurrentId, mReferenceBytes, mReferenceStart, mReferenceEnd - mReferenceStart);
      Diagnostic.developerLog("Moving reference buffer to: " + mCurrentName + "[" + (mReferenceStart + 1) + "," + (mReferenceStart + length) + ")");
      mCalibrator.setTemplate(mCurrentName, mReferenceStart, mReferenceBytes, length);
    }
  }

  /**
   * Return the calibrator used by this populator.
   * @return calibrator
   */
  public Calibrator calibrator() {
    return mCalibrator;
  }
}

