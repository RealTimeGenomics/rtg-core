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
    try (VcfReader r = VcfReader.openVcfReader(input, TabixIndexer.indexFileName(input).exists() ? region : null)) {

      // We have names available
      final List<String> sampleNames = r.getHeader().getSampleNames();
      if (r.getHeader().getNumberOfSamples() == 0) {
        throw new NoTalkbackSlimException("VCF file: " + input + " does not contain any samples.");
      } else if (r.getHeader().getNumberOfSamples() == 1) { // Assume the user asked for it for a reason, always just add it
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
