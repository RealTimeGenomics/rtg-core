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
package com.rtg.sam;

import com.rtg.util.machine.MachineType;

import net.sf.samtools.SAMReadGroupRecord;
import net.sf.samtools.SAMRecord;

/**
 */
public final class ReadGroupUtils {

  /** Read group attribute tag */
  public static final String RG_ATTRIBUTE = "RG";
  /** Name assigned to stats when no read group has been specified */
  public static final String UNKNOWN_RG = "unspecified";

  private ReadGroupUtils() { }

  /**
   * Get read group from <code>SAMRecord</code>
   * Will return 'unspecified' when read group is not set in record
   * @param record the SAM record to get the read group from
   * @return the read group
   */
  public static String getReadGroup(SAMRecord record) {
    final String rgAtt = record.getStringAttribute(RG_ATTRIBUTE);
    final String rgId;
    if (rgAtt == null) {
      rgId = UNKNOWN_RG;
    } else {
      rgId = "" + rgAtt;
    }
    return rgId;
  }

  /**
   * Get read group from <code>SamBamRecord</code>
   * Will return 'unspecified' when read group is not set in record
   * @param record the SAM record to get the read group from
   * @return the read group
   */
  public static String getReadGroup(SamBamRecord record) {
    final String rg = (String) record.getAttributeValue(RG_ATTRIBUTE);
    if (rg == null) {
      return UNKNOWN_RG;
    }
    return rg;
  }

  /**
   * Acquire machine type from read group if possible. Assumes platform has been set in read group
   * @param srgr the read group record.
   * @param paired if reads are paired or not
   * @return machine type if recognized, otherwise null
   */
  public static MachineType platformToMachineType(SAMReadGroupRecord srgr, boolean paired) {
    if (MachineType.COMPLETE_GENOMICS.compatiblePlatform(srgr.getPlatform())) {
      return MachineType.COMPLETE_GENOMICS;
    } else if (MachineType.FOURFIVEFOUR_PE.compatiblePlatform(srgr.getPlatform()) || MachineType.FOURFIVEFOUR_SE.compatiblePlatform(srgr.getPlatform())) {
      if (paired) {
        return MachineType.FOURFIVEFOUR_PE;
      } else {
        return MachineType.FOURFIVEFOUR_SE;
      }
    } else if (MachineType.ILLUMINA_PE.compatiblePlatform(srgr.getPlatform()) || MachineType.ILLUMINA_SE.compatiblePlatform(srgr.getPlatform())) {
      if (paired) {
        return MachineType.ILLUMINA_PE;
      } else {
        return MachineType.ILLUMINA_SE;
      }
    } else if (MachineType.IONTORRENT.compatiblePlatform(srgr.getPlatform())) {
      return MachineType.IONTORRENT;
    }
    return null;
  }
}
