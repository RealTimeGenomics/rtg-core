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

package com.rtg.simulation;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.reeltwo.jumble.annotations.TestClass;
import com.rtg.util.Pair;
import com.rtg.util.intervals.RegionRestriction;
import com.rtg.vcf.VcfReader;
import com.rtg.vcf.VcfRecord;
import com.rtg.vcf.VcfUtils;

/**
 * A class to wrap the offsets for a given mutated sample.
 */
@TestClass("com.rtg.simulation.MutatedReferenceReadNameParserTest")
public final class MutatedSampleOffsets {

  private final Map<String, Pair<MutatedOffsets, MutatedOffsets>> mOffsets = new HashMap<>();

  private MutatedSampleOffsets() { }

  /**
   * Create an offsets wrapper for the given sample.
   * @param phasedMutations the VCF file containing the phased mutations being compensated for.
   * @param sample the name of the sample being processed from the VCF file.
   * @throws IOException if a problem occurs when reading the VCF file.
   * @return the offsets for the sample.
   */
  public static MutatedSampleOffsets getOffsets(File phasedMutations, String sample) throws IOException {
    return getOffsets(phasedMutations, null, sample)[0];
  }

  /**
   * Create a set of offsets wrappers for the given samples.
   * @param phasedMutations the VCF file containing the phased mutations being compensated for.
   * @param region the region of the VCF file to construct a parser for.
   * @param samples the names of the samples to be processed from the VCF file. If empty, use all samples.
   * @throws IOException if a problem occurs when reading the VCF file.
   * @return the the offsets for the samples.
   */
  public static MutatedSampleOffsets[] getOffsets(File phasedMutations, RegionRestriction region, String... samples) throws IOException {
    try (VcfReader reader = VcfReader.openVcfReader(phasedMutations, region)) {
      final List<String> sampleNames = reader.getHeader().getSampleNames();
      final Map<String, Integer> sampleIndexMap = new HashMap<>(sampleNames.size());
      for (int i = 0; i < sampleNames.size(); ++i) {
        sampleIndexMap.put(sampleNames.get(i), i);
      }
      final MutatedSampleOffsets[] offsets = new MutatedSampleOffsets[samples.length];
      final int[] sampleIndexes = new int[samples.length];
      for (int i = 0; i < offsets.length; ++i) {
        offsets[i] = new MutatedSampleOffsets();
        sampleIndexes[i] = sampleIndexMap.get(samples[i]);
      }
      while (reader.hasNext()) {
        final VcfRecord rec = reader.next();
        final String ref = rec.getRefCall();
        final List<String> alleles = rec.getAltCalls();
        boolean indel = false;
        for (String al : alleles) {
          if (al.length() != ref.length()) {
            indel = true;
          }
        }
        if (!indel) {
          continue;
        }
        final List<String> gts = rec.getFormat(VcfUtils.FORMAT_GENOTYPE);
        for (int i = 0; i < samples.length; ++i) {
          final String gt = gts.get(sampleIndexes[i]);
          if (VcfUtils.isMissingGt(gt)) {
            continue;
          }
          final MutatedSampleOffsets mutatedReferenceOffsets = offsets[i];
          if (mutatedReferenceOffsets != null) {
            try {
              final int[] alleleIndexes = VcfUtils.splitGt(gt);
              if (alleleIndexes[0] > 0) {
                final int size = ref.length() - alleles.get(alleleIndexes[0] - 1).length();
                if (size != 0) {
                  final MutatedOffsets chromosomeOffsets = mutatedReferenceOffsets.getOrCreateOffsetPair(rec.getSequenceName()).getA();
                  if (size > 0) {
                    chromosomeOffsets.addDeletion(rec.getOneBasedStart(), size);
                  } else {
                    chromosomeOffsets.addInsertion(rec.getOneBasedStart(), -size);
                  }
                }
              }
              if (alleleIndexes.length == 2 && alleleIndexes[1] > 0) {
                final int size = ref.length() - alleles.get(alleleIndexes[1] - 1).length();
                if (size != 0) {
                  final MutatedOffsets chromosomeOffsets = mutatedReferenceOffsets.getOrCreateOffsetPair(rec.getSequenceName()).getB();
                  if (size > 0) {
                    chromosomeOffsets.addDeletion(rec.getOneBasedStart(), size);
                  } else {
                    chromosomeOffsets.addInsertion(rec.getOneBasedStart(), -size);
                  }
                }
              }
            } catch (IllegalArgumentException e) {
              //Overlapping variants, disable offsets
              offsets[i] = null;
            }
          }
        }
      }
      for (MutatedSampleOffsets sampleOffsets : offsets) {
        if (sampleOffsets != null) {
          for (Pair<MutatedOffsets, MutatedOffsets> chromosomeOffsets : sampleOffsets.mOffsets.values()) {
            chromosomeOffsets.getA().freeze();
            chromosomeOffsets.getB().freeze();
          }
        }
      }
      return offsets;
    }
  }

  private Pair<MutatedOffsets, MutatedOffsets> getOrCreateOffsetPair(String chromosomeName) {
    return mOffsets.computeIfAbsent(chromosomeName, k -> new Pair<>(new MutatedOffsets(), new MutatedOffsets()));
  }

  /**
   * Get the pair of offset structures associated with the given chromosome.
   * @param chromosome the chromosome to get the offsets for.
   * @return the pair of offsets for the given chromosome.
   */
  public Pair<MutatedOffsets, MutatedOffsets> get(String chromosome) {
    return mOffsets.get(chromosome);
  }

}
