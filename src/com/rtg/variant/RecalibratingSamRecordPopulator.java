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
package com.rtg.variant;

import java.io.IOException;
import java.util.Map;

import com.rtg.calibrate.Calibrator;
import com.rtg.reader.ReaderUtils;
import com.rtg.reader.SequencesReader;
import com.rtg.util.diagnostic.NoTalkbackSlimException;

import htsjdk.samtools.SAMRecord;

/**
 * Populator that passes through SAM records and performs calibration.
 */
public class RecalibratingSamRecordPopulator extends SamRecordPopulator {

  private static final byte[] EMPTY = new byte[0];
  private final Calibrator mCalibrator;
  private final SequencesReader mTemplate;
  private final byte[] mTemplateBytes;
  private final Map<String, Long> mSequenceNameMap;
  private String mLastName = "";

  /**
   * Populator for SAM records using an initializer. Can be given null to
   * skip initialization.
   * @param calibrator calibrator to use
   * @param template template reader
   */
  public RecalibratingSamRecordPopulator(final Calibrator calibrator, final SequencesReader template) {
    super();
    mCalibrator = calibrator;
    mTemplate = template;
    mTemplateBytes = new byte[(int) mTemplate.maxLength()];
    try {
      mSequenceNameMap = ReaderUtils.getSequenceNameMap(mTemplate);
    } catch (final IOException e) {
      // I don't think this should happen?
      throw new NoTalkbackSlimException("Problem getting sequence names from template");
    }
  }

  @Override
  public SAMRecord populate(final SAMRecord rec) {
    final String name = rec.getReferenceName();
    if (!name.equals(mLastName)) {
      if (!"*".equals(name)) {
        final Long seqId = mSequenceNameMap.get(name);
        if (seqId == null) {
          throw new NoTalkbackSlimException("Sequence " + name + " not found in template");  //user must have edited the sam file and cocked this up.
        }
        try {
          mTemplate.seek(seqId);
          final int length = mTemplate.readCurrent(mTemplateBytes);
          mCalibrator.setTemplate(mTemplateBytes, length);
        } catch (final IOException e) {
          throw new NoTalkbackSlimException("Failed to read sequence " + name + " from template");
        }
      } else {
        mCalibrator.setTemplate(EMPTY, 0);
      }
      mLastName = name;
    }
    mCalibrator.processRead(rec);
    return super.populate(rec);
  }

  /**
   * Return the calibrator used by this populator.
   * @return calibrator
   */
  public Calibrator calibrator() {
    return mCalibrator;
  }
}

