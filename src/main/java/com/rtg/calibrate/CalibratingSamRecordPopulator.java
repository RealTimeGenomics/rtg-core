/*
 * Copyright (c) 2018. Real Time Genomics Limited.
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the
 *    distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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

