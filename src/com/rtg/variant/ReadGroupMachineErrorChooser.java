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

import java.io.File;
import java.io.IOException;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;

import com.rtg.util.InvalidParamsException;
import com.rtg.util.StringUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.NoTalkbackSlimException;
import com.rtg.variant.util.VariantUtils;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMReadGroupRecord;

/**
 * Select appropriate machine errors for a given SAM record from its read group
 */
public class ReadGroupMachineErrorChooser implements MachineErrorChooserInterface {

  private final Map<String, MachineErrorParams> mMachineErrors = new HashMap<>();

  /**
   * Constructor
   * @param header header containing read group information
   * @throws IOException if an IO error occurs
   */
  public ReadGroupMachineErrorChooser(SAMFileHeader header) throws IOException {
    try {
      final List<SAMReadGroupRecord> groups = header.getReadGroups();
      if (groups.isEmpty()) {
        throw new NoTalkbackSlimException("No read groups found. Unable to determine machine error rate.");
      }
      addReadGroups(groups);
    } catch (final InvalidParamsException e) {
      throw new NoTalkbackSlimException(e.getMessage());
    }
  }

  /**
   * @param samInputs SAM files (with read groups).
   * @throws IOException when reading SAM headers.
   */
  public ReadGroupMachineErrorChooser(final File[] samInputs) throws IOException {
    try {
      for (final File file : samInputs) {
        try (SAMFileReader reader = new SAMFileReader(file)) {
          addReadGroups(reader.getFileHeader().getReadGroups());
        }
      }
    } catch (final InvalidParamsException e) {
      throw new NoTalkbackSlimException(e.getMessage());
    }
  }

  /**
   * @param samInputs SAM files (with read groups).
   * @throws IOException when reading SAM headers.
   */
  public ReadGroupMachineErrorChooser(final Collection<File> samInputs) throws IOException {
    this(samInputs.toArray(new File[samInputs.size()]));
  }

  private void addReadGroups(final List<SAMReadGroupRecord> groups) throws InvalidParamsException, IOException {
    if (groups.size() == 0) {
      throw new NoTalkbackSlimException("No read groups found. Unable to determine machine error rate.");
    }
    for (final SAMReadGroupRecord record : groups) {
      final String fPlatform = record.getPlatform();
      if (fPlatform != null) {
        if (fPlatform.equalsIgnoreCase("LS454")) {
          if (record.getPredictedMedianInsertSize() != null) {
            mMachineErrors.put(record.getId(), MachineErrorParams.builder("ls454_pe").create());
          } else {
            mMachineErrors.put(record.getId(), MachineErrorParams.builder("ls454_se").create());
          }
        } else {
          mMachineErrors.put(record.getId(), MachineErrorParams.builder(fPlatform.toLowerCase(Locale.ROOT)).create());
        }
      } else {
        throw new NoTalkbackSlimException("Read group: " + record.getId() + " has no specified platform. Unable to determine machine error rate.");
      }
      Diagnostic.developerLog("Machine errors for read group: " + record.getId() + StringUtils.LS + VariantUtils.dumpMachineErrors(mMachineErrors.get(record.getId())));
    }

  }

  @Override
  public MachineErrorParams machineErrors(VariantAlignmentRecord r) {
    final SAMReadGroupRecord rgr = r.getReadGroup();
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
}
