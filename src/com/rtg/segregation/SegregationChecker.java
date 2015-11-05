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

package com.rtg.segregation;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;

import com.reeltwo.jumble.annotations.TestClass;
import com.rtg.reference.ReferenceSequence;
import com.rtg.reference.Sex;
import com.rtg.util.ChiSquared;
import com.rtg.util.MathUtils;
import com.rtg.util.MultiSet;
import com.rtg.util.Pair;
import com.rtg.util.intervals.RangeList;
import com.rtg.util.Utils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.variant.format.VcfFormatField;
import com.rtg.vcf.VcfReader;
import com.rtg.vcf.VcfRecord;
import com.rtg.vcf.VcfUtils;
import com.rtg.vcf.VcfWriter;
import com.rtg.vcf.header.MetaType;
import com.rtg.vcf.header.SampleField;
import com.rtg.vcf.header.VcfHeader;
import com.rtg.vcf.header.VcfNumber;

/**
 */
@TestClass("com.rtg.segregation.SegregationCheckerCliTest")
public class SegregationChecker {

  private static final String PHASING_COMPATIBLE = "PHC";
  private static final String PHASING_REPAIRED = "PHR";
  private static final String PHASING_INCOMPATIBLE = "PHI";
  private static final String PHASING_OUTSIDE = "PHO";
  private static final String PHASING_INCOMPATIBLE_COUNT = "PHIC";
  private static final String PHASING_QUALITY = "PHQ";
  private static final String FILTER_PLOIDY_MISMATCH = "PME";
  private static final String ALT_FATHER_GT_FIX = "AFGT";
  private static final String ALT_MOTHER_GT_FIX = "AMGT";
  private static final String PHASE_GROUP_COUNTS = "PGC";
  private static final String PHASE_GROUP_DEPTH_RATIOS = "PGDR";
  private static final String PHASE_SET = "PS";

  private final VcfReader mReader;
  private final VcfWriter mWriter;
  private final Map<String, RangeList<PatternHolder>> mPatterns;
  private final Map<String, LinkedSet<PatternHolder>> mLinkedPatterns;
  private final boolean mRepairSimple;

  private final Map<Pair<Sex, String>, ReferenceSequence> mPloidyMap;
  private final int mSampleFather;
  private final int mSampleMother;
  private final Sex[] mSampleSex;
  private final List<String> mSampleNames;

  private PatternHolder mLastPatternHolder = null;
  private Integer mCurrPhaseGroup = null;
  private String mLastSeq = null;

  SegregationChecker(String father, String mother, VcfReader reader, VcfWriter writer, Map<String, RangeList<PatternHolder>> patterns, Map<Pair<Sex, String>, ReferenceSequence> ploidyMap, boolean repairSimple) {
    mReader = reader;
    mWriter = writer;
    mPatterns = patterns;
    mLinkedPatterns = new HashMap<>();
    for (final Entry<String, RangeList<PatternHolder>> el : mPatterns.entrySet()) {
      final LinkedSet<PatternHolder> ls = new LinkedSet<>();
      for (final RangeList.RangeData<PatternHolder> range : el.getValue().getRangeList()) {
        for (final PatternHolder pah : range.getMeta()) {
          ls.add(pah);
        }
      }
      mLinkedPatterns.put(el.getKey(), ls);
    }
    mPloidyMap = ploidyMap;
    mRepairSimple = repairSimple;
    final VcfHeader header = reader.getHeader();
    final Map<String, Integer> sampleIndexMap = new HashMap<>();
    mSampleNames = header.getSampleNames();
    for (int i = 0; i < mSampleNames.size(); i++) {
      sampleIndexMap.put(mSampleNames.get(i), i);
    }
    mSampleFather = sampleIndexMap.get(father);
    mSampleMother = sampleIndexMap.get(mother);
    final Map<String, Sex> sampleSexMap = new HashMap<>();
    for (final SampleField line : header.getSampleLines()) {
      sampleSexMap.put(line.getId(), line.getSex());
    }
    mSampleSex = SegregationVcfSearch.getSampleSex(mSampleNames, sampleSexMap, mSampleFather, mSampleMother);
  }

  void run() throws IOException {
    final Map<String, Integer> mismatchPloidyCounts = new LinkedHashMap<>();
    while (mReader.hasNext()) {
      final VcfRecord rec = mReader.next();
      try {
        processRecord(rec);
      } catch (final MismatchingPloidyException e) {
        Integer val = mismatchPloidyCounts.get(rec.getSequenceName());
        if (val == null) {
          val = 0;
        }
        mismatchPloidyCounts.put(rec.getSequenceName(), val + 1);
        rec.addFilter(FILTER_PLOIDY_MISMATCH);
      }
      mWriter.write(rec);
    }
    for (final Entry<String, Integer> seq : mismatchPloidyCounts.entrySet()) {
      Diagnostic.warning("There were " + seq.getValue() + " variants with genotypes that did not match the expected ploidy in chromosome " + seq.getKey());
    }
  }

  private FamilyGt familyFromRecord(VcfRecord rec) throws MismatchingPloidyException {
    assert mSampleFather != mSampleMother;
    final GType[] children = new GType[rec.getNumberOfSamples() - 2];
    final List<String> gts = rec.getFormatAndSample().get(VcfUtils.FORMAT_GENOTYPE);
    final GType father = new GType(gts.get(mSampleFather), mPloidyMap.get(new Pair<>(mSampleSex[0], rec.getSequenceName())).effectivePloidy(rec.getStart()));
    final GType mother = new GType(gts.get(mSampleMother), mPloidyMap.get(new Pair<>(mSampleSex[1], rec.getSequenceName())).effectivePloidy(rec.getStart()));
    int index = 0;
    for (int i = 0; i < gts.size(); i++) {
      if (i == mSampleFather || i == mSampleMother) {
        continue;
      }
      children[index] = new GType(gts.get(i), mPloidyMap.get(new Pair<>(mSampleSex[index + 2], rec.getSequenceName())).effectivePloidy(rec.getStart()));
      index++;
    }
    return new FamilyGt(rec.getSequenceName(), rec.getStart(), father, mother, children, false);
  }

  PatternHolder overRide(VcfRecord rec, String fatherGtOverride, String motherGtOverride, int childIndex, String childGtOverride) throws MismatchingPloidyException {
    assert mSampleFather != mSampleMother;
    final GType[] children = new GType[rec.getNumberOfSamples() - 2];
    final List<String> gts = rec.getFormatAndSample().get(VcfUtils.FORMAT_GENOTYPE);
    final GType father = new GType(fatherGtOverride == null ? gts.get(mSampleFather) : fatherGtOverride, mPloidyMap.get(new Pair<>(mSampleSex[0], rec.getSequenceName())).effectivePloidy(rec.getStart()));
    final GType mother = new GType(motherGtOverride == null ? gts.get(mSampleMother) : motherGtOverride, mPloidyMap.get(new Pair<>(mSampleSex[1], rec.getSequenceName())).effectivePloidy(rec.getStart()));
    int index = 0;
    for (int i = 0; i < gts.size(); i++) {
      if (i == mSampleFather || i == mSampleMother) {
        continue;
      }
      children[index] = new GType(childGtOverride == null || childIndex != index ? gts.get(i) : childGtOverride, mPloidyMap.get(new Pair<>(mSampleSex[index + 2], rec.getSequenceName())).effectivePloidy(rec.getStart()));
      index++;
    }
    final FamilyGt family = new FamilyGt(rec.getSequenceName(), rec.getStart(), father, mother, children, false);
    return new PatternHolder(family.pattern(), false);
  }

  void processRecord(final VcfRecord rec) throws IOException, MismatchingPloidyException {
    final String sequence = rec.getSequenceName();
    final List<PatternHolder> confident = getConfident(rec, sequence);
    if (confident == null) {
      rec.addFilter(PHASING_OUTSIDE);
      return;
    }
    if (!sequence.equals(mLastSeq)) {
      mLastPatternHolder = null;
      mCurrPhaseGroup = null;
      mLastSeq = sequence;
    }
    assert !((mLastPatternHolder == null) ^ (mCurrPhaseGroup == null));
    assert mLastSeq.equals(sequence);
    final VcfRecord repairedRecord = VcfReader.vcfLineToRecord(rec.toString());
    final PatternHolder recPhasing = new PatternHolder(familyFromRecord(rec).pattern(), false);
    final Set<String> possibleGts = getPossibleGts(rec.getAltCalls().size());
    PatternHolder selectedPatternChildren = null;
    PatternHolder selectedPatternParents = null;
    int minIncompatible = Integer.MAX_VALUE;
    String[] selectedIncompatible = null;
    for (final PatternHolder pattern : confident) {
      if (pattern.compatible(recPhasing)) {
        selectedPatternChildren = pattern;
        minIncompatible = 0;
        break;
      } else {
        for (final String alt : possibleGts) {
          try {
            final PatternHolder father = overRide(rec, alt, null, -1, null);
            final PatternHolder mother = overRide(rec, null, alt, -1, null);
            if (pattern.compatible(father) || pattern.compatible(mother)) {
              selectedPatternParents = pattern;
            }
          } catch (final MismatchingPloidyException e) {
            //TODO: perhaps set up alternates appropriate for differing ploidies?
          }
        }
        final String[] incompat = pattern.incompatibleChildren(recPhasing);
        int countIncompat = 0;
        for (final String compat : incompat) {
          if ("I".equals(compat)) {
            countIncompat++;
          }
        }
        if (countIncompat < minIncompatible) {
          selectedPatternChildren = pattern;
          selectedIncompatible = incompat;
          minIncompatible = countIncompat;
        }
      }
    }
    final int[] groupCounts = new int[4];
    final double[] groupDprSums = new double[4];
    final double[] groupDprAvgs = new double[4];
    int child = 0;
    for (int i = 0; i < rec.getNumberOfSamples(); i++) {
      if (mSampleFather == i || mSampleMother == i) {
        continue;
      }
      final int childGroup = selectedPatternChildren.group(child);
      groupCounts[childGroup]++;
      groupDprSums[childGroup] += rec.getFormatAndSample().get(VcfFormatField.DPR.name()) == null ? 0.0 : Double.valueOf(rec.getFormatAndSample().get(VcfFormatField.DPR.name()).get(i));
      child++;
    }
    for (int i = 0; i < groupCounts.length; i++) {
      if (groupCounts[i] > 0) {
        groupDprAvgs[i] = groupDprSums[i] / groupCounts[i];
      }
    }
    rec.addInfo(PHASE_GROUP_COUNTS, "" + groupCounts[0], "" + groupCounts[1], "" + groupCounts[2], "" + groupCounts[3]);
    rec.addInfo(PHASE_GROUP_DEPTH_RATIOS, Utils.realFormat(groupDprAvgs[0], 3), Utils.realFormat(groupDprAvgs[1], 3), Utils.realFormat(groupDprAvgs[2], 3), Utils.realFormat(groupDprAvgs[3], 3));
    final List<String> fatherAlts = new ArrayList<>();
    final List<String> motherAlts = new ArrayList<>();
    if (selectedPatternParents != null) {
      for (final String alt : possibleGts) {
        try {
          final PatternHolder father = overRide(rec, alt, null, -1, null);
          if (selectedPatternParents.compatible(father)) {
            fatherAlts.add(alt);
          }
        } catch (final MismatchingPloidyException e) {
          //TODO: perhaps set up alternates appropriate for differing ploidies?
        }
        try {
          final PatternHolder mother = overRide(rec, null, alt, -1, null);
          if (selectedPatternParents.compatible(mother)) {
            motherAlts.add(alt);
          }
        } catch (final MismatchingPloidyException e) {
          //TODO: perhaps set up alternates appropriate for differing ploidies?
        }
      }
      if (fatherAlts.size() > 0) {
        rec.addInfo(ALT_FATHER_GT_FIX, fatherAlts.toArray(new String[fatherAlts.size()]));
      }
      if (motherAlts.size() > 0) {
        rec.addInfo(ALT_MOTHER_GT_FIX, motherAlts.toArray(new String[motherAlts.size()]));
      }
    }
    if (selectedPatternChildren.compatible(recPhasing)) {
      rec.addInfo(PHASING_COMPATIBLE);
      phase(rec, selectedPatternChildren);
      calcPhq(rec);
    } else {
      rec.addFilter(PHASING_INCOMPATIBLE);
      rec.addInfo(PHASING_INCOMPATIBLE, selectedIncompatible);
      rec.addInfo(PHASING_INCOMPATIBLE_COUNT, Integer.toString(minIncompatible));
      if (mRepairSimple) {
        if (minIncompatible == 1) {
          int childIndex = -1;
          for (final String compat : selectedIncompatible) {
            childIndex++;
            if ("I".equals(compat)) {
              break;
            }
          }
          final List<String> childAlts = new ArrayList<>();
          for (final String alt : possibleGts) {
            try {
              final PatternHolder childAlt = overRide(repairedRecord, null, null, childIndex, alt);
              if (selectedPatternChildren.compatible(childAlt)) {
                childAlts.add(alt);
              }
            } catch (final MismatchingPloidyException e) {
              //TODO: perhaps set up alternates appropriate for differing ploidies?
            }
          }
          if (childAlts.size() == 1) {
            int sampleIndex = -1;
            for (int i = 0; i < mSampleNames.size(); i++) {
              sampleIndex++;
              if (i == mSampleFather || i == mSampleMother) {
                continue;
              } else if (childIndex == 0) {
                break;
              }
              childIndex--;
            }
            repairRecord(rec, repairedRecord, sampleIndex, childAlts.get(0), selectedPatternChildren);
          }
        } else if (fatherAlts.size() == 1 && motherAlts.size() == 0) {
          repairRecord(rec, repairedRecord, mSampleFather, fatherAlts.get(0), selectedPatternParents);
        } else if (motherAlts.size() == 1 && fatherAlts.size() == 0) {
          repairRecord(rec, repairedRecord, mSampleMother, motherAlts.get(0), selectedPatternParents);
        }
      }
    }
  }

  List<PatternHolder> getConfident(final VcfRecord rec, final String sequence) {
    final RangeList<PatternHolder> ranges = mPatterns.get(sequence);
    if (ranges == null) {
      return null;
    }
    return ranges.find(rec.getStart());
  }

  static Set<String> getPossibleGts(final int altSize) {
    final Set<String> possibleGts = new HashSet<>();
    for (int i = 0; i < altSize; i++) {
      for (int j = i; j < altSize; j++) {
        possibleGts.add(i + "/" + j);
      }
    }
    return possibleGts;
  }

  static VcfHeader modifyHeader(final VcfHeader header, boolean repairSimple) {
    header.addInfoField(PHASING_COMPATIBLE, MetaType.FLAG, new VcfNumber("0"), "The phasing of children in this variant is compatible with known phasing pattern");
    header.addInfoField(PHASING_INCOMPATIBLE, MetaType.CHARACTER, VcfNumber.DOT, "The phasing of children in this variant is incompatible with known phasing pattern, 'C' -> consistent or 'I' -> inconsistent for each child in the order of children in the samples.");
    header.addInfoField(PHASING_INCOMPATIBLE_COUNT, MetaType.INTEGER, VcfNumber.ONE, "Count of the minimum number of inconsistent children.");
    header.addInfoField(PHASING_QUALITY, MetaType.INTEGER, VcfNumber.ONE, "Phred-scaled probability that the phasing consistency would have been obtained by chance.");
    header.addInfoField(ALT_FATHER_GT_FIX, MetaType.STRING, VcfNumber.DOT, "List of alternate GT values for father that would make call consistent.");
    header.addInfoField(ALT_MOTHER_GT_FIX, MetaType.STRING, VcfNumber.DOT, "List of alternate GT values for mother that would make call consistent.");
    header.addInfoField(PHASE_GROUP_COUNTS, MetaType.INTEGER, new VcfNumber("4"), "The number of children in each phasing group (00, 01, 10, 11)");
    header.addInfoField(PHASE_GROUP_DEPTH_RATIOS, MetaType.FLOAT, new VcfNumber("4"), "The average DPR of children in each phasing group (00, 01, 10, 11)");
    header.addFilterField(PHASING_INCOMPATIBLE, "This variant has a phasing incompatibility");
    header.addFilterField(PHASING_OUTSIDE, "This variant was outside the regions of known phasing");
    header.addFilterField(FILTER_PLOIDY_MISMATCH, "This variant should be ignored as the genotype ploidy of some of the samples did not match the expected ploidy");
    header.addFormatField(PHASE_SET, MetaType.INTEGER, VcfNumber.ONE, "Phase set for the genotype");

    if (repairSimple) {
      header.addFilterField(PHASING_REPAIRED, "This variant should be ignored as it has been replaced by a variant with a repair for phasing incompatibility");
      header.addInfoField(PHASING_REPAIRED, MetaType.STRING, new VcfNumber("2"), "The phasing of children in this variant has been repaired by altering one samples GT, the two strings are the sample that was repaired and the old GT");
    }
    return header;
  }

  void repairRecord(VcfRecord original, VcfRecord rec, int sampleIndex, String newGt, PatternHolder curr) throws IOException, MismatchingPloidyException {
    final String sampleName = mSampleNames.get(sampleIndex);
    final String oldGt = rec.getFormatAndSample().get(VcfUtils.FORMAT_GENOTYPE).get(sampleIndex);
    rec.setFormatAndSample(VcfUtils.FORMAT_GENOTYPE, newGt, sampleIndex);
    rec.addInfo(PHASING_REPAIRED, sampleName, oldGt);
    rec.addInfo(PHASING_COMPATIBLE);
    original.addFilter(PHASING_REPAIRED);
    phase(rec, curr);
    calcPhq(rec);
    mWriter.write(rec);
  }

  private void calcPhq(VcfRecord rec) {
    if (!rec.isFiltered()) {
      final MultiSet<String> codes = new MultiSet<>();
      for (String gtStr : rec.getFormatAndSample().get(VcfUtils.FORMAT_GENOTYPE)) {
        final int[] gtint = VcfUtils.splitGt(gtStr);
        final String code;
        if (gtint.length == 1) {
          code = "" + gtint[0];
        } else if (gtint.length == 2 && gtint[1] < gtint[0]) {
          code = "" + gtint[1] + "|" + gtint[0];
        } else {
          code = "" + gtint[0] + "|" + gtint[1];
        }
        codes.add(code);
      }
      rec.addInfo(PHASING_QUALITY, "" + (int) calcPhq(codes));
    }
  }

  /**
   * Implement calculation of the probability that d different genotypes,
   * occurring across both parents and children with counts given as arguments,
   * will agree with phasing by chance.
   *
   * @param counts contain counts for each of the distinct genotypes
   * @return the phred-scaled probability
   */
  static double calcPhq(MultiSet<String> counts) {
    final int d = counts.keySet().size();
    final double est = 1.0 / d;
    final long nsum = counts.totalCount();
    double lprob = 0;
    for (final String key : counts.keySet()) {
      lprob += ChiSquared.lgamma(counts.get(key) + est);
    }
    lprob -= ChiSquared.lgamma(nsum + 1);
    lprob -= d * ChiSquared.lgamma(est);
    return MathUtils.lnToPhred(lprob);
  }

  /**
   * Add phasing to rec given the pattern in curr.
   * @param rec the, possibly copied, record.
   * @param curr the pattern to use for phasing.
   * @throws MismatchingPloidyException shouldnt happen.
   */
  private void phase(VcfRecord rec, PatternHolder curr) throws MismatchingPloidyException {
    //System.err.println("phase");
    final FamilyGt family = familyFromRecord(rec);
    final PatternArray pattern = curr.pattern();
    final Labelling labels = new Labelling(family, pattern);
    final GType father = family.father();
    final GType mother = family.mother();
    if (father.isSingleAllele() && mother.isSingleAllele()) {
      return;
    }
    for (int i = 0, j = -1; i < mSampleNames.size(); i++) {
      if (i == mSampleFather) {
        if (father.isSingleAllele()) {
          continue;
        }
        phaseSet(rec, i, labels.phaseFather(), curr);
      } else if (i == mSampleMother) {
        if (mother.isSingleAllele()) {
          continue;
        }
        phaseSet(rec, i, labels.phaseMother(), curr);
      } else {
        //child
        j++;
        if (family.child(j).isSingleAllele()) {
          continue;
        }
        phaseSet(rec, i, labels.phaseChild(j), curr);
      }
    }
    rec.padFormatAndSample(PHASE_SET);
  }

  private void phaseSet(final VcfRecord rec, final int i, final GType phase, final PatternHolder curr) {
    if (phase == null) {
      return;
    }
    if (i == mSampleFather || i == mSampleMother) {
      if (mLastPatternHolder == null || isNew(curr)) {
        mCurrPhaseGroup = rec.getOneBasedStart();
        mLastPatternHolder = curr;
      }
      rec.setFormatAndSample(PHASE_SET, mCurrPhaseGroup.toString(), i);
    }
    final String phStr = phase.a() + "|" + phase.b();
    //System.err.println(phStr);
    rec.setFormatAndSample(VcfUtils.FORMAT_GENOTYPE, phStr, i);
  }

  private boolean isNew(final PatternHolder curr) {
    final Iterator<PatternHolder> it = mLinkedPatterns.get(mLastSeq).iterator(mLastPatternHolder, curr);
    while (it.hasNext()) {
      if (it.next().isNew()) {
        return true;
      }
    }
    return false;
  }

}
