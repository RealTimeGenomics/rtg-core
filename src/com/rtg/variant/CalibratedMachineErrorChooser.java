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
import java.util.concurrent.ConcurrentHashMap;

import com.rtg.calibrate.Calibrator;
import com.rtg.sam.ReadGroupUtils;
import com.rtg.util.InvalidParamsException;
import com.rtg.util.StringUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.NoTalkbackSlimException;
import com.rtg.util.machine.MachineType;
import com.rtg.variant.util.VariantUtils;

import htsjdk.samtools.SAMReadGroupRecord;

/**
 * Machine error chooser that uses calibration files.
 */
public class CalibratedMachineErrorChooser implements MachineErrorChooserInterface {

  private static final boolean CG_BYPASS_HACK = true; //Boolean.valueOf(System.getProperty("cg-calibration-bypass", "true"));

  private static AbstractMachineErrorParams sCompleteParams = null;

  private static synchronized AbstractMachineErrorParams getDefaultCompleteParams(MachineType mt) {
    if (sCompleteParams == null) {
      try {
        sCompleteParams = MachineErrorParams.builder(mt.priors()).create();
      } catch (final IOException e) {
        throw new RuntimeException("Could not load built-in complete genomics error rates");
      } catch (final InvalidParamsException e) {
        // Converted here because this is a build-problem and we don't want to add InvalidParamsException to callers of this constructor
        throw new RuntimeException("Bad built-in complete genomics error rates");
      }
    }
    return sCompleteParams;
  }


  private final Calibrator mCalibrator;

  private final Map<String, AbstractMachineErrorParams> mReadGroupMachineErrorParams;

  /**
   * Constructor
   * @param c calibrator to choose on
   */
  public CalibratedMachineErrorChooser(Calibrator c) {
    mCalibrator = c;
    mReadGroupMachineErrorParams = new ConcurrentHashMap<>();
  }

  @Override
  public AbstractMachineErrorParams machineErrors(VariantAlignmentRecord r) {
    final SAMReadGroupRecord rg = r.getReadGroup();
    if (rg == null) {
      throw new NoTalkbackSlimException("Read group required in SAM file");
    }
    final String rgId = rg.getId();
    AbstractMachineErrorParams cal = mReadGroupMachineErrorParams.get(rgId);
    if (cal == null) {
      final MachineType mt = ReadGroupUtils.platformToMachineType(rg, r.isReadPaired());
      if (mt == MachineType.COMPLETE_GENOMICS && CG_BYPASS_HACK) {
        Diagnostic.developerLog("CG calibration bypass enabled, using default CG errors");
        cal = getDefaultCompleteParams(mt);
      } else {
        cal = new CalibratedMachineErrorParams(mt, mCalibrator, rgId);
      }
      Diagnostic.developerLog("Machine errors for read group: " + rgId + StringUtils.LS + VariantUtils.dumpMachineErrors(cal));
      mReadGroupMachineErrorParams.put(rgId, cal);
    }
    return cal;
  }

}
