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

import java.lang.reflect.Constructor;
import java.lang.reflect.InvocationTargetException;
import java.util.ArrayList;
import java.util.List;

import com.rtg.launcher.globals.CoreGlobalFlags;
import com.rtg.launcher.globals.GlobalFlags;
import com.rtg.mode.DNA;
import com.rtg.mode.DnaUtils;
import com.rtg.util.SeparateClassLoader;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.SlimException;
import com.rtg.util.intervals.SequenceNameLocusSimple;
import com.rtg.variant.bayes.Hypotheses;
import com.rtg.variant.realign.AbstractAllPaths;
import com.rtg.variant.realign.AlignmentEnvironment;
import com.rtg.variant.realign.AlignmentEnvironmentGenome;
import com.rtg.variant.realign.AlignmentEnvironmentGenomeSubstitution;
import com.rtg.variant.realign.AllPaths;
import com.rtg.variant.realign.Environment;
import com.rtg.variant.realign.EnvironmentCombined;
import com.rtg.variant.realign.RealignParams;
import com.rtg.variant.realign.RealignParamsGenome;
import com.rtg.variant.realign.ScoreFastUnderflow;
import com.rtg.variant.realign.ScoreMatrix;
import com.rtg.variant.util.arithmetic.PossibilityArithmetic;

/**
 */
public class ComplexTemplate extends SequenceNameLocusSimple {

  /*
   * The following block of code is being used to allow two separate class definitions
   * for classes that the java JIT can optimise quite heavily in another part of the code
   * if not used directly in this code. This allows it to maintain the speed it used to have
   * instead of becoming approximately 40% slower due to the optimisations removed to allow
   * both locations to use the same code.
   *
   * TODO: THIS WAS TRUE WITH JAVA 6, 7, and 8 -- EVALUATE IF STILL NEEDED
   */
  /*******************START BLOCK******************/
  private static final boolean HOTSPOT_HACK = GlobalFlags.getBooleanValue(CoreGlobalFlags.COMPLEX_CALLER_HOTSPOT_HACK);
  private static final Constructor<AllPaths> SCORE_FAST_CONSTRUCTOR;
  private static final Constructor<Environment> ENVIRONMENT_COMBINED_CONSTRUCTOR;
  static {
    Diagnostic.developerLog("Classloading hack = " + HOTSPOT_HACK);
    final SeparateClassLoader loader =  new SeparateClassLoader(ScoreFastUnderflow.class, ScoreMatrix.class, AbstractAllPaths.class, EnvironmentCombined.class);
    try {
      @SuppressWarnings("unchecked")
      final Class<AllPaths> sfClazz = (Class<AllPaths>) loader.loadClass(ScoreFastUnderflow.class.getName());
      SCORE_FAST_CONSTRUCTOR = sfClazz.getConstructor(RealignParams.class);
      @SuppressWarnings("unchecked")
      final Class<Environment> ecClazz = (Class<Environment>) loader.loadClass(EnvironmentCombined.class.getName());
      ENVIRONMENT_COMBINED_CONSTRUCTOR = ecClazz.getConstructor(AlignmentEnvironment.class, int.class, int.class, AlignmentEnvironment.class);
    } catch (final ClassNotFoundException | NoSuchMethodException | SecurityException e) {
      throw new SlimException(e);
    }
  }

  private static AllPaths initScoreFastUnderflow(RealignParams params) {
    if (HOTSPOT_HACK) {
      try {
        return SCORE_FAST_CONSTRUCTOR.newInstance(params);
      } catch (final IllegalArgumentException | InvocationTargetException | IllegalAccessException | InstantiationException e) {
        throw new SlimException(e);
      }
    } else {
      return new ScoreFastUnderflow(params);
    }
  }

  private static Environment initEnvironmentCombined(AlignmentEnvironment samEnv, int zeroBasedStartPos, int maxShift, AlignmentEnvironment templateEnv) {
    if (HOTSPOT_HACK) {
      try {
        return ENVIRONMENT_COMBINED_CONSTRUCTOR.newInstance(samEnv, zeroBasedStartPos, maxShift, templateEnv);
      } catch (final IllegalArgumentException | InvocationTargetException | IllegalAccessException | InstantiationException e) {
        throw new SlimException(e);
      }
    } else {
      return new EnvironmentCombined(samEnv, zeroBasedStartPos, maxShift, templateEnv);
    }
  }
  /*******************END BLOCK********************/


  private final byte[] mTemplate;
  private final byte[] mReplaceBytes;
  private final String mReplaceString;
  private DescriptionComplex mDescription;
  private PossibilityArithmetic mArithmetic;
  private int mRefHyp;
  private int mRefTransitionIndex;
  private double[][] mTransitionProbsLn;

  /**
   * @param template the bytes of the full template.
   * @param refName name of reference sequence.
   * @param start position in the template of the region to be replaced (0 based inclusive).
   * @param end  position in the template of the region to be replaced (0 based exclusive).
   */
  public ComplexTemplate(byte[] template, String refName, int start, int end) {
    super(refName, start, end);
    mTemplate = template;
    final int length = end - start;
    mReplaceBytes = new byte[length];
    final StringBuilder sb = new StringBuilder();
    for (int i = 0, j = start; i < length; ++i, ++j) {
      final byte b = mTemplate[j];
      mReplaceBytes[i] = b;
      sb.append(DnaUtils.getBase(b));
    }
    mReplaceString = sb.toString();
  }

  private static int findReferenceAllele(DescriptionComplex description, String reference) {
    for (int i1 = 0; i1 < description.size(); ++i1) {
      if (description.name(i1).equals(reference)) {
        return i1;
      }
    }
    return Hypotheses.NO_HYPOTHESIS;
  }

  /**
   * Get the underlying full template bytes.
   * @return the underlying template bytes.
   */
  public byte[] templateBytes() {
    return mTemplate;
  }

  /**
   * Get the nucleotides in the replacement region as bytes.
   * @return the bytes in the replacement region.
   */
  public byte[] replaceBytes() {
    return mReplaceBytes;
  }

  /**
   * Get the string of nucleotides corresponding to the region to be replaced.
   * @return the string of nucleotides corresponding to the region to be replaced ("" for zero length).
   */
  public String replaceString() {
    return mReplaceString;
  }

  /**
   * Set the information needed to compute all-paths priors
   * @param description the description to be used during complex calling
   * @param arithmetic the arithmetic
   */
  public void setComplexContext(DescriptionComplex description, PossibilityArithmetic arithmetic) {
    mDescription = description;
    mArithmetic = arithmetic;
    mRefHyp = findReferenceAllele(description, replaceString());
    computeTransitionMatrix();
  }

  private void computeTransitionMatrix() {
    final AllPaths sm = initScoreFastUnderflow(RealignParamsGenome.SINGLETON);
    final int hypExtension = Math.max(5, Math.max(description().maxLength() + 1, getLength() + 1));
    final int start = getStart() - hypExtension;
    final int end = getStart() + hypExtension;
    final int e0Hx = getEnd() - hypExtension;

    // Build up all alignment environments that we will need
    final int n = description().size();
    final List<AlignmentEnvironment> envs = new ArrayList<>(n); // All alleles to which we will perform alignment (description first,
    final List<Integer> lengths = new ArrayList<>(n);
    for (int i = 0; i < n; ++i) {
      final AlignmentEnvironment env;
      final int len;
      if (i == refHyp()) {
        env = new AlignmentEnvironmentGenome(start, end, templateBytes());
        len = getEnd() - getStart();
      } else {
        final byte[] hypDna = DNA.stringDNAtoByte(description().name(i));
        len = hypDna.length;
        env = new AlignmentEnvironmentGenomeSubstitution(start, end, this, hypDna);
      }
      envs.add(env);
      lengths.add(len);
    }

    final int refAllele;
    if (refHyp() == Hypotheses.NO_HYPOTHESIS) { // Reference wasn't included in the set of hypotheses (perhaps due to Ns)
      refAllele = envs.size();
      envs.add(new AlignmentEnvironmentGenome(start, end, templateBytes()));
      lengths.add(getLength());
    } else {
      refAllele = refHyp();
    }
    mRefTransitionIndex = refAllele;

    // Compute all-paths transition matrix
    final double[][] transitionProbsLn = new double[envs.size()][n];
    for (int i = 0; i < envs.size(); ++i) {
      final AlignmentEnvironment aei = envs.get(i);
      final int iLen = lengths.get(i);
      for (int j = 0; j < n; ++j) {
        final AlignmentEnvironment aej = envs.get(j);
        final int jLen = lengths.get(j);
        final int alignStart = e0Hx - (jLen + iLen) / 2;
        final int maxShift = (jLen + iLen + 1) / 2;
        final Environment env = initEnvironmentCombined(aej, alignStart, maxShift, aei);
        sm.setEnv(env);
        transitionProbsLn[i][j] = sm.totalScoreLn();
      }
    }
    mTransitionProbsLn = transitionProbsLn;
  }

  /**
   * @return the description to be used during complex calling
   */
  public DescriptionComplex description() {
    return mDescription;
  }

  /**
   * Gets the matrix of transition probabilities between each description name. There may be an additional row in the matrix
   * if the reference sequence is not an element of the description.
   * @return the matrix of transition probabilities between each description name.
   */
  public double[][] transitionProbsLn() {
    return mTransitionProbsLn;
  }

  PossibilityArithmetic arithmetic() {
    return mArithmetic;
  }

  int refHyp() {
    return mRefHyp;
  }

  int refTransitionIndex() {
    return mRefTransitionIndex;
  }

}
