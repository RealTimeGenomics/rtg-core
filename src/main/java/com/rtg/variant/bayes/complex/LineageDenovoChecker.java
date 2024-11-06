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

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import com.rtg.relation.LineageLookup;
import com.rtg.util.StringUtils;
import com.rtg.variant.Variant;
import com.rtg.variant.VariantSample;
import com.rtg.variant.util.VariantUtils;

/**
 * De novo checker that supports cell lineage de novo checking.
 */
public class LineageDenovoChecker implements DenovoChecker {

  final LineageLookup mLineage;

  /**
   * @param lineage array of length equal to the number of samples, giving the the
   * sample id of the predecessor (parent or originating sample) for each sample
   */
  public LineageDenovoChecker(LineageLookup lineage) {
    mLineage = lineage;
  }

  @Override
  public boolean isDenovo(Variant variant, int sample) {
    final VariantSample derived = variant.getSample(sample);
    final int originalId = mLineage.getOriginal(sample);
    final VariantSample original = variant.getSample(originalId);
    final List<String> childAlleles = new ArrayList<>(Arrays.asList(StringUtils.split(derived.getName(), VariantUtils.COLON)));
    for (final String s : StringUtils.split(original.getName(), VariantUtils.COLON)) {
      if (!childAlleles.remove(s)) {
        return true;
      }
    }
    return false;
  }
}
