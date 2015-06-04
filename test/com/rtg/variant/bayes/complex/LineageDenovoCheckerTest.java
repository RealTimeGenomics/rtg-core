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
package com.rtg.variant.bayes.complex;

import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import com.rtg.reference.Ploidy;
import com.rtg.relation.LineageLookup;
import com.rtg.util.StringUtils;
import com.rtg.variant.Variant;
import com.rtg.variant.VariantLocus;
import com.rtg.variant.VariantSample;
import com.rtg.variant.bayes.snp.DescriptionCommon;

import junit.framework.TestCase;

/**
 */
public class LineageDenovoCheckerTest extends TestCase {
  public static boolean isDenovo(DenovoChecker checker, String ref, String[] ... samples) {
    final VariantLocus locus = new VariantLocus("chr", 0, ref.length(), ref, 'N');
    final VariantSample[] variantSamples = new VariantSample[samples.length];
    final Set <String> alleles = new HashSet<>();
    for (String[] sample : samples) {
      if (sample != null) {
        final List<String> sampleList = Arrays.asList(sample);
        alleles.addAll(sampleList);
      }
    }
    alleles.add(ref);
    final DescriptionCommon description = TrimmingTest.getDescription(alleles.toArray(new String[alleles.size()]));
    for (int i = 0; i < samples.length; i++) {
      final List<String> sampleList = samples[i] == null ? null : Arrays.asList(samples[i]);
      variantSamples[i] = samples[i] == null ? null : TrimmingTest.getVariantSample(samples[i].length == 1 ? Ploidy.HAPLOID : Ploidy.DIPLOID, StringUtils.join(":", sampleList), false, 5.0, VariantSample.DeNovoStatus.IS_DE_NOVO, 0.0, description);

    }
    final Variant v = new Variant(locus, variantSamples);
    return checker.isDenovo(v, samples.length - 1);

  }
  public void testDenovoCorrect() {
    final LineageDenovoChecker checker = new LineageDenovoChecker(new LineageLookup(-1, 0));
    assertTrue(isDenovo(checker, "T", new String[] {"T", "T"}, new String[] {"T", "A"}));
    assertTrue(isDenovo(checker, "T", new String[] {"T", "T"}, new String[] {"T", "A"}));
    assertTrue(isDenovo(checker, "A", new String[] {"TT", "TT"}, new String[] {"TT", "A"}));
    assertTrue(isDenovo(checker, "A", new String[]{"T", "C"}, new String[]{"T", "T"}));

    assertFalse(isDenovo(checker, "T", new String[] {"T", "T"}, new String[] {"T", "T"}));
    assertFalse(isDenovo(checker, "T", new String[] {"T", "C"}, new String[] {"T", "C"}));
    assertFalse(isDenovo(checker, "T", new String[]{"T", "C"}, new String[]{"C", "T"}));
  }
}
