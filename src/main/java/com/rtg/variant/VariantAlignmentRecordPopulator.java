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

import java.util.HashMap;

import com.rtg.launcher.globals.CoreGlobalFlags;
import com.rtg.launcher.globals.GlobalFlags;
import com.rtg.sam.HomopolymerUtils;
import com.rtg.util.Populator;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.NoTalkbackSlimException;

import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;

/**
 * Populator.
 */
public class VariantAlignmentRecordPopulator implements Populator<VariantAlignmentRecord> {

  private final HashMap<String, Integer> mGenomeToInteger;
  private final MachineErrorChooserInterface mChooser;
  private final int mMinBaseQuality;
  private final boolean mMaskHomopolymer;
  private long mMaskedRecords;
  private long mTotalRecords;

  /**
   * Populator.
   * @param chooser Machine error chooser. Will be used to recalibrate qualities at record creation time
   * @param minBaseQuality minimum read base quality in phred scale
   * @param genomeNames list of genome names
   */
  public VariantAlignmentRecordPopulator(MachineErrorChooserInterface chooser, int minBaseQuality, final String... genomeNames) {
    mChooser = chooser;
    mGenomeToInteger = new HashMap<>(genomeNames.length);
    mMinBaseQuality = minBaseQuality;
    for (int i = 0; i < genomeNames.length; ++i) {
      final String s = genomeNames[i];
      mGenomeToInteger.put(s, i);
    }
    mMaskHomopolymer = GlobalFlags.getBooleanValue(CoreGlobalFlags.VARIANT_MASK_HOMOPOLYMER);
    mMaskedRecords = 0;
    mTotalRecords = 0;
  }

  @Override
  public VariantAlignmentRecord populate(final SAMRecord rec) {
    if (mMaskHomopolymer) {
      ++mTotalRecords;
      if (HomopolymerUtils.maskHomoPolymer(rec)) {
        ++mMaskedRecords;
      }
      if (mTotalRecords % 1000000 == 0) {
        Diagnostic.developerLog("Masked " + mMaskedRecords + " out of " + mTotalRecords + " (" + (100.0 * mMaskedRecords / mTotalRecords) + "%)");
      }
    }
    if (mGenomeToInteger.size() > 0) {
      final SAMReadGroupRecord readGroup = rec.getReadGroup();
      if (readGroup == null) {
        throw new NoTalkbackSlimException("Encountered a SAM record with no read group information: " + rec.getSAMString());
      }
      final Integer genome = mGenomeToInteger.get(readGroup.getSample());
      if (genome == null) {
        throw new NoTalkbackSlimException("Could not determine sample from SAM record (check read group information against expected samples): " + rec.getSAMString());
      }
      try {
        return new VariantAlignmentRecord(rec, genome, mChooser, mMinBaseQuality);
      } catch (IllegalArgumentException e) {
        return null;
      }
    } else {
      try {
        return new VariantAlignmentRecord(rec, 0, mChooser, mMinBaseQuality);
      } catch (IllegalArgumentException e) {
        return null;
      }
    }
  }

  @Override
  public VariantAlignmentRecord overflow(int position, int length) {
    return VariantAlignmentRecord.overflow(position, length);
  }
}

