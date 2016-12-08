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
      final Map<String, Integer> sampleIndexMap = new HashMap<>();
      final List<String> sampleNames = reader.getHeader().getSampleNames();
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
    Pair<MutatedOffsets, MutatedOffsets> ret = mOffsets.get(chromosomeName);
    if (ret == null) {
      ret = new Pair<>(new MutatedOffsets(), new MutatedOffsets());
      mOffsets.put(chromosomeName, ret);
    }
    return ret;
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
