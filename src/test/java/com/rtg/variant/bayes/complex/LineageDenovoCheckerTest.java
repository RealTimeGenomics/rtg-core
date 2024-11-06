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
import com.rtg.variant.util.VariantUtils;

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
    final DescriptionCommon description = SimpleDecomposerTest.getDescription(alleles.toArray(new String[0]));
    for (int i = 0; i < samples.length; ++i) {
      final List<String> sampleList = samples[i] == null ? null : Arrays.asList(samples[i]);
      variantSamples[i] = samples[i] == null ? null : SimpleDecomposerTest.getVariantSample(samples[i].length == 1 ? Ploidy.HAPLOID : Ploidy.DIPLOID, StringUtils.join("" + VariantUtils.COLON, sampleList), false, 5.0, VariantSample.DeNovoStatus.IS_DE_NOVO, 0.0, description);

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
