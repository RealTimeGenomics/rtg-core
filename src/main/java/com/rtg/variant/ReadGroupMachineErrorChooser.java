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
