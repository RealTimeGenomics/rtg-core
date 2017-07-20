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

package com.rtg.visualization;

import java.io.File;
import java.io.IOException;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import com.reeltwo.jumble.annotations.TestClass;
import com.rtg.mode.DnaUtils;
import com.rtg.reader.SequencesReader;
import com.rtg.sam.ReadGroupUtils;
import com.rtg.sam.SamFilterParams;
import com.rtg.sam.SamReadingContext;
import com.rtg.sam.SamRecordPopulator;
import com.rtg.sam.SamRegionRestriction;
import com.rtg.sam.SamUtils;
import com.rtg.sam.ThreadedMultifileIterator;
import com.rtg.util.IntegerOrPercentage;
import com.rtg.util.SingletonPopulatorFactory;
import com.rtg.util.diagnostic.NoTalkbackSlimException;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;

/**
 * Helper functions
 */
public final class SamHelper {

  private SamHelper() { }

  static final Comparator<SAMRecord> POSITION_COMP = new PositionComparator();
  static final Comparator<SAMRecord> RG_COMP = new ReadGroupsComparator();
  static final Comparator<SAMRecord> SAMPLE_COMP = new SamplesComparator();

  @TestClass("com.rtg.visualization.AviewModelTest")
  private static class PositionComparator implements Comparator<SAMRecord>, Serializable {
    @Override
    public int compare(SAMRecord s1, SAMRecord s2) {
      int res = s1.getAlignmentStart() - s2.getAlignmentStart();
      if (res == 0) {
        res = s1.getAlignmentEnd() - s2.getAlignmentEnd();
      }
      return res;
    }
  }

  @TestClass("com.rtg.visualization.AviewModelTest")
  private static class TwoStringComparator extends PositionComparator {

    protected int compare(SAMRecord s1, String str1, SAMRecord s2, String str2) {
      final int compare;
      if (str1 != null && str2 != null) {
        compare = str1.compareTo(str2);
      } else {
        if (str1 == null) {
          if (str2 == null) {
            compare = 0; //rg1 == null rg2 == null
          } else {
            compare = 1; //rg1 == null rg2 != null
          }
        } else {
          compare = -1; //rg1 != null rg2 == null;
        }
      }
      if (compare == 0) {
        return super.compare(s1, s2);
      }
      return compare;
    }
  }

  @TestClass("com.rtg.visualization.AviewModelTest")
  private static class ReadGroupsComparator extends TwoStringComparator {
    @Override
    public int compare(SAMRecord s1, SAMRecord s2) {
      return super.compare(s1, s1.getStringAttribute(ReadGroupUtils.RG_ATTRIBUTE),
        s2, s2.getStringAttribute(ReadGroupUtils.RG_ATTRIBUTE));
    }
  }
  @TestClass("com.rtg.visualization.AviewModelTest")
  private static class SamplesComparator extends TwoStringComparator {
    @Override
    public int compare(SAMRecord s1, SAMRecord s2) {
      return super.compare(s1, s1.getReadGroup() != null ? s1.getReadGroup().getSample() : null,
        s2, s2.getReadGroup() != null ? s2.getReadGroup().getSample() : null);
    }
  }

  static String getCorrectComplement(boolean forward, final String readString) {
    if (forward) {
      return readString;
    }
    return DnaUtils.reverseComplement(readString);
  }

  static ArrayList<SAMRecord> loadAlignments(final AviewParams params, SequencesReader reader) throws IOException {
    final ArrayList<SAMRecord> records = new ArrayList<>();
    final List<File> files = Arrays.asList(params.alignmentsFiles());
    if (files.size() == 0) {
      return records;
    }
    final SAMFileHeader header = SamUtils.getUberHeader(reader, files);
    if (header.getSequence(params.sequenceName()) == null) {
      throw new NoTalkbackSlimException("Unable to apply region \"" + params.sequenceName() + "\" as the specified sequence name does not exist in input files.");
    }

    // If they have asked to restrict by sample, create read group white list
    Set<String> validRGs = null;
    final String[] samples = params.wantedSamples();
    if (samples != null) {
      validRGs = new HashSet<>();
      for (final SAMReadGroupRecord rg : header.getReadGroups()) {
        final String sm = rg.getSample();
        for (final String sample : samples) {
          if (sample.equals(sm)) {
            validRGs.add(rg.getId());
          }
        }
      }
    }

    final SamFilterParams sfp = SamFilterParams.builder()
      .minMapQ(params.minMapQ())
      .maxMatedAlignmentScore(new IntegerOrPercentage(params.maxMatedAlignmentScore()))
      .maxUnmatedAlignmentScore(new IntegerOrPercentage(params.maxUnmatedAlignmentScore()))
      .restriction(new SamRegionRestriction(params.region())).create();
    final SingletonPopulatorFactory<SAMRecord> pf = new SingletonPopulatorFactory<>(new SamRecordPopulator());
    final SamReadingContext context = new SamReadingContext(files, 2, sfp, header, reader);
    try (final ThreadedMultifileIterator<SAMRecord> it = new ThreadedMultifileIterator<>(context, pf)) {
      while (it.hasNext()) {
        final SAMRecord r = it.next();
        if (validIhScore(r, params.maxIhScore())
            && validReadGroup(validRGs, r.getStringAttribute(ReadGroupUtils.RG_ATTRIBUTE))) {
          records.add(r);
        }
      }
    }
    return records;
  }

  static boolean validReadGroup(Set<String> validRGs, String rg) {
    return validRGs == null || rg == null || validRGs.contains(rg);
  }

  static boolean validIhScore(SAMRecord r, int maxScore) {
    final Integer nh = SamUtils.getNHOrIH(r);
    if (nh == null) {
      return true;
    }
    return nh < maxScore;
  }
}
