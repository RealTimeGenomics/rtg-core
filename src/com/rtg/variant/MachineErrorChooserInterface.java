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

import com.rtg.util.machine.MachineType;
import com.rtg.variant.realign.RealignParams;

import htsjdk.samtools.SAMReadGroupRecord;

/**
 */
public interface MachineErrorChooserInterface {

  /**
   * @param rg read group to select errors from
   * @param readPaired is the read paired
   * @return appropriate machine errors for read
   */
  PhredScaler machineErrors(SAMReadGroupRecord rg, boolean readPaired);

  /**
   * @param rg read group of interest
   * @param readPaired is the read paired
   * @return realign parameters reflecting the machine errors observed in the read group
   */
  RealignParams realignParams(SAMReadGroupRecord rg, boolean readPaired);

  /**
   * @param rg read group of interest
   * @param readPaired is the read paired
   * @return type of machine this read group is from
   */
  MachineType machineType(SAMReadGroupRecord rg, boolean readPaired);
}
