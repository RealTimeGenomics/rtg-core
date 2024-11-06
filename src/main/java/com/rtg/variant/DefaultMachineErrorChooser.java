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

import com.rtg.util.InvalidParamsException;
import com.rtg.util.StringUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.machine.MachineType;
import com.rtg.variant.realign.RealignParams;
import com.rtg.variant.realign.RealignParamsImplementation;
import com.rtg.variant.util.VariantUtils;

import htsjdk.samtools.SAMReadGroupRecord;

/**
 * Machine Error Chooser used when machine specified on the command line
 */
public class DefaultMachineErrorChooser implements MachineErrorChooserInterface {

  private final MachineErrorParams mMachineError;
  private final RealignParams mRealignParams;
  private String mPlatform = null;
  private boolean mWarned;

  /**
   * @throws InvalidParamsException whenever.
   */
  public DefaultMachineErrorChooser() {
    this(MachineErrorParams.builder().create());
  }

  /**
   * @param machineError default to be used for all SAM records.
   */
  public DefaultMachineErrorChooser(final MachineErrorParams machineError) {
    mMachineError = machineError;
    mRealignParams = new RealignParamsImplementation(mMachineError);
    Diagnostic.developerLog("Machine errors for all read groups: " + StringUtils.LS + VariantUtils.dumpMachineErrors(mMachineError));
  }

  /**
   * @param errorName name of errors to be used as default to be used for all SAM records.
   * @throws IOException when reading the errors file.
   */
  public DefaultMachineErrorChooser(final String errorName) throws IOException {
    mMachineError = MachineErrorParams.builder(errorName).create();
    mRealignParams = new RealignParamsImplementation(mMachineError);
    Diagnostic.developerLog("Machine errors for all read groups: " + errorName + StringUtils.LS + VariantUtils.dumpMachineErrors(mMachineError));
  }

  static final String PLATFORM_WARNING =
    "You are using a single set of machine error rates when combining reads from different platforms. Usually using platform specific error rates will give better results";

  @Override
  public MachineErrorParams machineErrors(SAMReadGroupRecord rg, boolean readPaired) {
    // Don't bother checking platforms if we've already warned
    if (!mWarned) {
      if (rg != null) {
        final String platform = rg.getPlatform();
        if (mPlatform == null) {
          mPlatform = platform;
        } else {
          if (!mPlatform.equals(platform)) {
            Diagnostic.warning(PLATFORM_WARNING);
            mWarned = true;
          }
        }
      }
    }
    return mMachineError;
  }

  @Override
  public RealignParams realignParams(SAMReadGroupRecord rg, boolean readPaired) {
    return mRealignParams;
  }

  @Override
  public MachineType machineType(SAMReadGroupRecord rg, boolean readPaired) {
    return machineErrors(rg, readPaired).machineType();
  }
}
