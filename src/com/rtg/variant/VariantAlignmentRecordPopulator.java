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

import java.util.HashMap;

import com.rtg.util.Populator;
import com.rtg.util.diagnostic.NoTalkbackSlimException;

import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;

/**
 * Populator.
 */
public class VariantAlignmentRecordPopulator implements Populator<VariantAlignmentRecord> {

  private final HashMap<String, Integer> mGenomeToInteger;

  /**
   * Populator.
   * @param genomeNames list of genome names
   */
  public VariantAlignmentRecordPopulator(final String... genomeNames) {
    mGenomeToInteger = new HashMap<>(genomeNames.length);
    for (int i = 0; i < genomeNames.length; i++) {
      final String s = genomeNames[i];
      mGenomeToInteger.put(s, i);
    }
  }

  @Override
  public VariantAlignmentRecord populate(final SAMRecord rec) {
    if (mGenomeToInteger.size() > 0) {
      final SAMReadGroupRecord readGroup = rec.getReadGroup();
      if (readGroup == null) {
        throw new NoTalkbackSlimException("Encountered a SAM record with no read group information: " + rec.getSAMString());
      }
      final Integer genome = mGenomeToInteger.get(readGroup.getSample());
      if (genome == null) {
        throw new NoTalkbackSlimException("Could not determine sample from SAM record (check read group information against expected samples): " + rec.getSAMString());
      }
      return new VariantAlignmentRecord(rec, genome);
    } else {
      return new VariantAlignmentRecord(rec, 0);
    }
  }

  @Override
  public VariantAlignmentRecord overflow(int position, int length) {
    return VariantAlignmentRecord.overflow(position, length);
  }
}

