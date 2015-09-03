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
import java.util.Collections;
import java.util.Comparator;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import com.reeltwo.jumble.annotations.TestClass;
import com.rtg.mode.DnaUtils;
import com.rtg.sam.ReadGroupUtils;
import com.rtg.sam.SamFilterParams;
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

  static String getCorrectComplement(boolean forward, final String readString) {
    if (forward) {
      return readString;
    }
    return DnaUtils.reverseComplement(readString);
  }

  static ArrayList<SAMRecord> loadAlignments(final AviewParams params) throws IOException {
    final ArrayList<SAMRecord> records = new ArrayList<>();
    final List<File> files = Arrays.asList(params.alignmentsFiles());
    if (files.size() == 0) {
      return records;
    }
    final SAMFileHeader header = SamUtils.getUberHeader(files);
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

    final SingletonPopulatorFactory<SAMRecord> pf = new SingletonPopulatorFactory<>(new SamRecordPopulator());
    try (final ThreadedMultifileIterator<SAMRecord> it = new ThreadedMultifileIterator<>(files, 2, pf,
      SamFilterParams.builder()
        .minMapQ(params.minMapq())
        .maxMatedAlignmentScore(new IntegerOrPercentage(params.maxMatedAlignmentScore()))
        .maxUnmatedAlignmentScore(new IntegerOrPercentage(params.maxUnmatedAlignmentScore()))
        .restriction(new SamRegionRestriction(params.region())).create(), header)) {
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

  static void sortAlignments(ArrayList<SAMRecord> records) {
    Collections.sort(records, new Comparator<SAMRecord>() {
        @Override
        public int compare(SAMRecord s1, SAMRecord s2) {
          int res = s1.getAlignmentStart() - s2.getAlignmentStart();
          if (res == 0) {
            res = s1.getAlignmentEnd() - s2.getAlignmentEnd();
          }
          return res;
        }
      });
  }

  @TestClass("com.rtg.visualization.AviewModelTest")
  private static class AlignmentsWithReadGroupsComparator implements Comparator<SAMRecord>, Serializable {
    @Override
    public int compare(SAMRecord s1, SAMRecord s2) {
      final String readGroup1 = s1.getStringAttribute(ReadGroupUtils.RG_ATTRIBUTE);
      final String readGroup2 = s2.getStringAttribute(ReadGroupUtils.RG_ATTRIBUTE);
      final int compare;
      if (readGroup1 != null && readGroup2 != null) {
        compare = readGroup1.compareTo(readGroup2);
      } else {
        if (readGroup1 == null) {
          if (readGroup2 == null) {
            compare = 0; //rg1 == null rg2 == null
          } else {
            compare = 1; //rg1 == null rg2 != null
          }
        } else {
          compare = -1; //rg1 != null rg2 == null;
        }
      }
      if (compare == 0) {
        int res = s1.getAlignmentStart() - s2.getAlignmentStart();
        if (res == 0) {
          res = s1.getAlignmentEnd() - s2.getAlignmentEnd();
        }
        return res;
      }
      return compare;
    }
  }

  static void sortAlignmentsWithReadGroups(ArrayList<SAMRecord> records) {
    Collections.sort(records, new AlignmentsWithReadGroupsComparator());
  }

  static boolean validIhScore(SAMRecord r, int maxScore) {
    final Integer nh = SamUtils.getNHOrIH(r);
    if (nh == null) {
      return true;
    }
    return nh < maxScore;
  }
}
