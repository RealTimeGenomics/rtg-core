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

import com.rtg.util.InvalidParamsException;
import com.rtg.util.StringUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.NoTalkbackSlimException;
import com.rtg.variant.util.VariantUtils;

import net.sf.samtools.SAMReadGroupRecord;

/**
 * Machine Error Chooser used when machine specified on the command line
 */
public class DefaultMachineErrorChooser implements MachineErrorChooserInterface {

  private final MachineErrorParams mMachineError;
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
    Diagnostic.developerLog("Machine errors for all read groups: " + StringUtils.LS + VariantUtils.dumpMachineErrors(mMachineError));
  }

  /**
   * @param errorName name of errors to be used as default to be used for all SAM records.
   * @throws IOException when reading the errors file.
   */
  public DefaultMachineErrorChooser(final String errorName) throws IOException {
    try {
      mMachineError = MachineErrorParams.builder(errorName).create();
      Diagnostic.developerLog("Machine errors for all read groups: " + errorName + StringUtils.LS + VariantUtils.dumpMachineErrors(mMachineError));
    } catch (final InvalidParamsException e) {
      throw new NoTalkbackSlimException(e.getMessage());
    }
  }

  static final String PLATFORM_WARNING =
    "You are using a single set of machine error rates when combining reads from different platforms. Usually using platform specific error rates will give better results";

  @Override
  public MachineErrorParams machineErrors(final VariantAlignmentRecord r) {
    // Don't bother checking platforms if we've already warned
    if (!mWarned) {
      final SAMReadGroupRecord rg = r.getReadGroup();
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
}
