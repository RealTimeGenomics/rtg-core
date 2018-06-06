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
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import com.rtg.tabix.TabixIndexer;
import com.rtg.util.intervals.RegionRestriction;
import com.rtg.util.diagnostic.NoTalkbackSlimException;
import com.rtg.vcf.VcfReader;
import com.rtg.vcf.VcfRecord;
import com.rtg.vcf.VcfUtils;

/**
 */
final class VariantHelper {

  private VariantHelper() { }

  /**
   * @param snpmap destination map for putting loaded SNPs
   * @param input the variant file
   * @param region specify the region for which SNPs will be loaded
   * @param wantedSampleNames only include SNPs corresponding to these samples
   * @throws IOException if there is a problem reading the variants
   */
  static void loadSnpRange(Map<String, ArrayList<AviewVariant>> snpmap, File input, RegionRestriction region, String... wantedSampleNames) throws IOException {
    // First look up the array of sample indexes we want to get from this file
    final ArrayList<Integer> actualSamples = new ArrayList<>();
    final ArrayList<String> actualSampleNames = new ArrayList<>();
    final File index = new File(input + TabixIndexer.TABIX_EXTENSION);
    try (VcfReader r = VcfReader.openVcfReader(input, index.exists() ? region : null)) {

      // We have names available
      final List<String> sampleNames = r.getHeader().getSampleNames();
      if (r.getHeader().getNumberOfSamples() == 1) { // Assume the user asked for it for a reason, always just add it
        actualSamples.add(0);
        final String sampleName = sampleNames.get(0);
        actualSampleNames.add(sampleName);
        snpmap.put(sampleName, new ArrayList<>());
      } else {
        // We have multiple samples, if the user didn't ask for a sample, print the list of availables and exit
        if (wantedSampleNames == null || wantedSampleNames.length == 0) {
          final StringBuilder sb = new StringBuilder("No sample name specified but multiple samples available. Please select from the samples available:");
          for (final String sampleName : sampleNames) {
            sb.append(" ").append(sampleName);
          }
          throw new NoTalkbackSlimException(sb.toString());
        }
        for (final String wanted : wantedSampleNames) {
          for (int i = 0; i < sampleNames.size(); ++i) {
            if (wanted.equals(sampleNames.get(i))) {
              actualSamples.add(i);
              actualSampleNames.add(wanted);
              snpmap.put(wanted, new ArrayList<>());
              break;
            }
          }
        }
      }

      while (r.hasNext()) {
        final VcfRecord current = r.next();
        //System.err.println("Got variant at " + (current.getPosition()) + ". range (" + start + "," + end + "]");
        // Explicitly check the position in case this is from an untabixed source.
        if (!region.overlaps(current)) {
          continue;
        }
        for (int i = 0; i < actualSamples.size(); ++i) {
          final String sampleName = actualSampleNames.get(i);
          final int sampleIndex = VcfUtils.getSampleIndexOrDie(r.getHeader(), sampleName, input.getPath());
          // Skip non-calls
          if (!VcfUtils.hasDefinedVariantGt(current, sampleIndex)) {
            continue;
          }
          final AviewVariant v = new AviewVariant(current, sampleIndex);
          final ArrayList<AviewVariant> snp = snpmap.get(sampleName);
          snp.add(v);
        }

      }
    }
  }

}
