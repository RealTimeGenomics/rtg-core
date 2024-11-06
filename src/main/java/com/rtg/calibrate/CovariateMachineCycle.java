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

  static int getReadLength(SAMRecord rec) {
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
