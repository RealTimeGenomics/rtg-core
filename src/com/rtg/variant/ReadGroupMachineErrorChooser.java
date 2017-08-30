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
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.rtg.sam.ReadGroupUtils;
import com.rtg.util.StringUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.NoTalkbackSlimException;
import com.rtg.util.machine.MachineType;
import com.rtg.variant.realign.RealignParams;
import com.rtg.variant.realign.RealignParamsImplementation;
import com.rtg.variant.util.VariantUtils;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;

/**
 * Select appropriate machine errors for a given SAM record from its read group
 */
public class ReadGroupMachineErrorChooser implements MachineErrorChooserInterface {

  private final Map<String, MachineErrorParams> mMachineErrors = new HashMap<>();
  private final Map<String, RealignParams> mRealignParams = new HashMap<>();

  /**
   * Constructor
   * @param header header containing read group information
   * @throws IOException if an IO error occurs
   */
  public ReadGroupMachineErrorChooser(SAMFileHeader header) throws IOException {
    addReadGroups(header.getReadGroups());
  }

  private void addReadGroups(final List<SAMReadGroupRecord> groups) throws IOException {
    if (groups.isEmpty()) {
      throw new NoTalkbackSlimException("No read groups found. Unable to determine machine error rate. Try explicitly specifying machine type");
    }
    for (final SAMReadGroupRecord record : groups) {
      final String fPlatform = record.getPlatform();
      if (fPlatform != null) {
        final MachineType mt = ReadGroupUtils.platformToMachineType(record, record.getPredictedMedianInsertSize() != null);
        if (mt != null) {
          final MachineErrorParams me = MachineErrorParams.builder(mt.priors()).create();
          mMachineErrors.put(record.getId(), me);
          mRealignParams.put(record.getId(), new RealignParamsImplementation(me));
        } else {
          throw new NoTalkbackSlimException("Read group: " + record.getId() + " has unrecognized platform. Unable to determine machine error rate. Try explicitly specifying machine type");
        }
      } else {
        throw new NoTalkbackSlimException("Read group: " + record.getId() + " has no specified platform. Unable to determine machine error rate. Try explicitly specifying machine type");
      }
      Diagnostic.developerLog("Machine errors for read group: " + record.getId() + StringUtils.LS + VariantUtils.dumpMachineErrors(mMachineErrors.get(record.getId())));
    }
  }

  @Override
  public MachineErrorParams machineErrors(SAMReadGroupRecord rgr, boolean readPaired) {
    final String rg = rgr == null ? null : rgr.getId();
    if (rg == null) {
      throw new NoTalkbackSlimException("Sam record had no read group attribute, but header read groups were supplied.");
    }

    final MachineErrorParams me = mMachineErrors.get(rg);
    if (me == null) {
      throw new NoTalkbackSlimException("Sam record referenced read group \"" + rg + "\" which was not found in the header.");
    }
    return me;
  }

  @Override
  public RealignParams realignParams(SAMReadGroupRecord rgr, boolean readPaired) {
    final String rg = rgr == null ? null : rgr.getId();
    if (rg == null) {
      throw new NoTalkbackSlimException("Sam record had no read group attribute, but header read groups were supplied.");
    }

    final RealignParams rp = mRealignParams.get(rg);
    if (rp == null) {
      throw new NoTalkbackSlimException("Sam record referenced read group \"" + rg + "\" which was not found in the header.");
    }
    return rp;
  }

  @Override
  public MachineType machineType(SAMReadGroupRecord rg, boolean readPaired) {
    return machineErrors(rg, readPaired).machineType();
  }
}
