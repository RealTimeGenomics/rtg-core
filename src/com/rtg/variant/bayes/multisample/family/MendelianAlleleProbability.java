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

package com.rtg.variant.bayes.multisample.family;

import com.rtg.variant.bayes.Code;

/**
 */
public interface MendelianAlleleProbability {

  /**
   * @param code description of encoding used for nucleotides.
   * @param father the hypothesis from the father's sample.
   * @param mother the hypothesis from the mother's sample. Can be -1 when there is no relevant code (eg mothers contribution to Y chromosome).
   * @param child the hypothesis from a child's sample.
   * @return the probability <code>M</code>.  See the <code>multiScoring.tex</code> documentation for Consistent Mendelian pedigrees.
   */
  double probabilityLn(Code code, int father, int mother, int child);

  /**
   * @param code description of encoding used for nucleotides.
   * @param father the hypothesis from the father's sample.
   * @param mother the hypothesis from the mother's sample. Can be -1 when there is no relevant code (eg mothers contribution to Y chromosome).
   * @param child the hypothesis from a child's sample.
   * @return true if the hypotheses are no mendelian consistent.
   */
  boolean isDenovo(Code code, int father, int mother, int child);
}
