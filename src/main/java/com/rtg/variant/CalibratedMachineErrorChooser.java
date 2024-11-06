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
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;

import com.rtg.calibrate.Calibrator;
import com.rtg.sam.ReadGroupUtils;
import com.rtg.util.InvalidParamsException;
import com.rtg.util.Pair;
import com.rtg.util.StringUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.NoTalkbackSlimException;
import com.rtg.util.machine.MachineType;
import com.rtg.variant.realign.RealignParams;
import com.rtg.variant.realign.RealignParamsImplementation;
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
        throw new RuntimeException("Could not load built-in complete genomics error rates", e);
      } catch (final InvalidParamsException e) {
        // Converted here because this is a build-problem and we don't want to add InvalidParamsException to callers of this constructor
        throw new RuntimeException("Bad built-in complete genomics error rates", e);
      }
    }
    return sCompleteParams;
  }


  private final Calibrator mCalibrator;

  private final Map<String, Pair<AbstractMachineErrorParams, RealignParams>> mReadGroupMachineErrorParams;

  /**
   * Constructor
   * @param c calibrator to choose on
   */
  public CalibratedMachineErrorChooser(Calibrator c) {
    mCalibrator = c;
    mReadGroupMachineErrorParams = new ConcurrentHashMap<>();
  }

  private Pair<AbstractMachineErrorParams, RealignParams> lookup(SAMReadGroupRecord rg, boolean readPaired) {
    if (rg == null) {
      throw new NoTalkbackSlimException("Read group required in SAM file");
    }
    final String rgId = rg.getId();
    Pair<AbstractMachineErrorParams, RealignParams> cr = mReadGroupMachineErrorParams.get(rgId);
    if (cr == null) {
      MachineType mt = ReadGroupUtils.platformToMachineType(rg, readPaired);
      if (mt == null) {
        Diagnostic.warning("Read group " + rg.getId() + " does not contain a recognized platform, assuming generic");
        mt = MachineType.GENERIC;
      }
      final AbstractMachineErrorParams cal;
      if (mt == MachineType.COMPLETE_GENOMICS && CG_BYPASS_HACK) {
        Diagnostic.developerLog("CG calibration bypass enabled, using default CG errors");
        cal = getDefaultCompleteParams(mt);
      } else {
        cal = new CalibratedMachineErrorParams(mt, mCalibrator, rgId);
      }
      Diagnostic.developerLog("Machine errors for read group: " + rgId + StringUtils.LS + VariantUtils.dumpMachineErrors(cal));
      cr = new Pair<>(cal, new RealignParamsImplementation(cal));
      mReadGroupMachineErrorParams.put(rgId, cr);
    }
    return cr;
  }

  @Override
  public AbstractMachineErrorParams machineErrors(SAMReadGroupRecord rg, boolean readPaired) {
    return lookup(rg, readPaired).getA();
  }

  @Override
  public RealignParams realignParams(SAMReadGroupRecord rg, boolean readPaired) {
    return lookup(rg, readPaired).getB();
  }

  @Override
  public MachineType machineType(SAMReadGroupRecord rg, boolean readPaired) {
    return machineErrors(rg, readPaired).machineType();
  }
}
