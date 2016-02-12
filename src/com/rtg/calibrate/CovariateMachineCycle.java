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

import java.util.HashMap;
import java.util.Locale;
import java.util.Map;

import com.rtg.sam.ReadGroupUtils;
import com.rtg.util.machine.MachineOrientation;
import com.rtg.util.machine.MachineType;

import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;

/**
 * Covariate variable for the position along the read as sequenced.  Will automatically
 * grow if longer reads are encountered.
 */
public final class CovariateMachineCycle extends CovariateImpl {

  /**
   * @param maxReadLen expected maximum length of all reads.
   */
  public CovariateMachineCycle(int maxReadLen) {
    super(CovariateEnum.MACHINECYCLE.name().toLowerCase(Locale.ROOT), maxReadLen);
  }

  @Override
  public String name() {
    return super.name() + ":" + size();
  }

  private final Map<SAMReadGroupRecord, MachineOrientation> mReadGroupToOrientation = new HashMap<>();

  private MachineOrientation getMachineOrientation(final SAMRecord rec) {
    final SAMReadGroupRecord rg = rec.getReadGroup();
    final MachineOrientation mo = mReadGroupToOrientation.get(rg);
    if (mo != null) {
      return mo;
    }
    final MachineType machineType = rg == null ? null : ReadGroupUtils.platformToMachineType(rg, rec.getReadPairedFlag());
    MachineOrientation mo2 = machineType == null ? null : machineType.orientation();
    if (mo2 == null) {
      mo2 = MachineOrientation.FR; // Assume an Illumina-style machine
    }
    mReadGroupToOrientation.put(rg, mo2);
    return mo2;
  }

  private boolean getOrientation(final SAMRecord rec) {
    final boolean rc = rec.getReadNegativeStrandFlag();
    switch (getMachineOrientation(rec)) {
      case ANY: // For purposes of machine cycle assume ANY is actually FR
      case FR:
        return rc;
      case RF:
        return !rc;
      case TANDEM:
        return false;
      default:
        throw new InternalError();
    }
  }

  @Override
  public int value(SAMRecord sam, CalibratorCigarParser parser) {
    final int length = sam.getReadLength();
    if (length >= newSize()) {
      setNewSize(length + 1);
    }
    final int readPos = parser.getReadPosition();
    return getOrientation(sam) ? length - readPos - 1 : readPos;
  }

  @Override
  public CovariateEnum getType() {
    return CovariateEnum.MACHINECYCLE;
  }
}
