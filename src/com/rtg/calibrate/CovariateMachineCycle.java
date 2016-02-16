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

import java.util.Locale;

import com.rtg.reader.CgUtils;
import com.rtg.sam.ReadGroupUtils;
import com.rtg.util.machine.MachineType;

import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;

/**
 * Covariate variable for the position along the read as sequenced.  Will automatically
 * grow as longer and longer reads are encountered.
 */
public final class CovariateMachineCycle extends CovariateImpl {

  /**
   * @param maxReadLen expected maximum length of all reads (can be 0).
   */
  public CovariateMachineCycle(int maxReadLen) {
    super(CovariateEnum.MACHINECYCLE.name().toLowerCase(Locale.ROOT), maxReadLen);
  }

  @Override
  public String name() {
    return super.name() + ":" + size();
  }

  private int getReadLength(SAMRecord rec) {
    // Complicated by Complete Genomics reads where the rec.getReadLength() is not the true flattened read length
    if (rec.getHeader() != null) {
      final SAMReadGroupRecord rg = rec.getReadGroup();
      if (rg != null) {
        final MachineType machineType = ReadGroupUtils.platformToMachineType(rg, rec.getReadPairedFlag());
        if (machineType == MachineType.COMPLETE_GENOMICS_2) {
          return CgUtils.CG2_RAW_READ_LENGTH;
        } else if (machineType == MachineType.COMPLETE_GENOMICS) {
          return CgUtils.CG_RAW_READ_LENGTH;
        }
      }
    }
    return rec.getReadLength();
  }

  @Override
  public int value(SAMRecord sam, CalibratorCigarParser parser) {
    final int length = getReadLength(sam);
    final int readPos = parser.getReadPosition();
    if (readPos >= newSize()) {
      // Try to avoid doing too many resizes, by stepping directly out to the read length.
      // Unfortunately, getReadLength is not the flattened length of a CG read, so the
      // actual read length can exceeds the length reported by the SAM record.
      setNewSize(Math.max(length, readPos + 1));
    }
    final boolean rc = sam.getReadNegativeStrandFlag();
    final int machineCycle = rc ? length - readPos - 1 : readPos;
    assert machineCycle >= 0 && machineCycle < newSize() : "pos=" + readPos + " len=" + length + " " + sam.getSAMString();
    return machineCycle;
  }

  @Override
  public CovariateEnum getType() {
    return CovariateEnum.MACHINECYCLE;
  }
}
