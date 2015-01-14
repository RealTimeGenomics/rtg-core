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

package com.rtg.variant.eval;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.Callable;

import com.rtg.util.Pair;
import com.rtg.util.intervals.RegionRestriction;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.vcf.VcfReader;
import com.rtg.vcf.VcfRecord;
import com.rtg.vcf.VcfUtils;
import com.rtg.vcf.header.VcfHeader;

/**
 * A callable which loads VcfRecords for given template
 */
public class VcfRecordTabixCallable implements Callable<LoadedVariants> {

  private final File mInput;
  private final String mSampleName;
  private final String mTemplateName;
  private final int mTemplateLength;
  private final VariantSetType mType;
  private final RocSortValueExtractor mExtractor;
  private final boolean mPassOnly;
  private final boolean mSquashPloidy;
  private final int mMaxLength;

  VcfRecordTabixCallable(File file, String sample, Pair<String, Integer> templateNameLength, VariantSetType type, RocSortValueExtractor extractor, boolean passOnly, boolean squashPloidy, int maxLength) {
    mInput = file;
    mSampleName = sample;
    mTemplateName = templateNameLength.getA();
    mTemplateLength = templateNameLength.getB();
    mType = type;
    mExtractor = extractor;
    mPassOnly = passOnly;
    mSquashPloidy = squashPloidy;
    mMaxLength = maxLength;
  }

  @Override
  public LoadedVariants call() throws Exception {
    int skipped = 0;
    final List<DetectedVariant> list = new ArrayList<>();
    try (VcfReader reader = VcfReader.openVcfReader(mInput, new RegionRestriction(mTemplateName))) {
      final VcfHeader header = reader.getHeader();
      final String label = mType == VariantSetType.BASELINE ? "baseline" : "calls";
      final int sampleId = VcfUtils.getSampleIndexOrDie(header, mSampleName, label);
      DetectedVariant last = null;
      while (reader.hasNext()) {
        final VcfRecord rec = reader.next();

        // Skip non-variant, SV, and possibly fail variants
        if (VcfUtils.skipRecordForSample(rec, sampleId, mPassOnly)) {
          continue;
        }

        // Skip variants that are too long (these cause problems during evaluation)
        int length = rec.getRefCall().length();
        for (String alt : rec.getAltCalls()) {
          length = Math.max(alt.length(), length);
        }
        if (mMaxLength > -1 && length > mMaxLength) {
          Diagnostic.userLog("Variant in " + label + " at " + rec.getSequenceName() + ":" + rec.getOneBasedStart() + " exceeds maximum length, skipping.");
          skipped++;
          continue;
        }
        
        // Skip variants with starts falling outside the expected length of the template sequence
        if (mTemplateLength >= 0 && rec.getStart() >= mTemplateLength) {
          Diagnostic.userLog("Variant in " + label + " at " + rec.getSequenceName() + ":" + rec.getOneBasedStart() + " starts outside the length of the reference sequence (" + mTemplateLength + ").");
          skipped++;
          continue;
        }

        // Skip overlapping variants (this detection code was moved here from Path)
        final DetectedVariant v = new DetectedVariant(rec, sampleId, mExtractor, mSquashPloidy);
        if (last != null) {
          // There is tricky case where an
          // insertion event occurs adjacent to a snp event. There might
          // still be other situations not properly covered here.
          final boolean atAnInsert = (v.getStart() == v.getEnd()) ^ (last.getStart() == last.getEnd());
          if ((atAnInsert && last.getStart() > v.getStart()) || (!atAnInsert && (last.getStart() >= v.getStart() || last.getEnd() > v.getStart()))) {
            Diagnostic.userLog("Overlapping variants aren't supported, skipping current variant from " + label + ".\nPrevious variant: " + last + "\nCurrent variant:  " + v);
            skipped++;
            continue;
          }
        }
        last = v;
        list.add(v);
      }
    }
    return new LoadedVariants(list, skipped);
  }
}

