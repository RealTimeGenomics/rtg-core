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

package com.rtg.variant.bayes.multisample;

import java.io.IOException;
import java.util.List;

import com.rtg.launcher.globals.GlobalFlags;
import com.rtg.launcher.globals.CoreGlobalFlags;
import com.rtg.reference.SexMemo;
import com.rtg.relation.Family;
import com.rtg.variant.VariantParams;
import com.rtg.variant.bayes.Code;
import com.rtg.variant.bayes.ModelInterface;

/**
 */
public final class Utils {

  /**
   * Minimum required coverage reference=N locations in order to produce a call. With minimum of 1,
   * cases of large regions tend to trigger a complex call with only a single (probably wrong) hypothesis. Increasing
   * this reduces the size of the complex region to contain more alternative hypotheses.
   */
  public static final int MIN_DEPTH_FOR_N_CALL = GlobalFlags.getIntegerValue(CoreGlobalFlags.CALLER_N_MIN_DEPTH);

  private Utils() { }

  private static void mark(final boolean[] mark, final Code code, final int hyp) {
    if (hyp != -1) {
      final int a = code.a(hyp);
      mark[a] = true;
      if (!code.homozygous(hyp)) {
        final int b = code.b(hyp);
        mark[b] = true;
      }
    }
  }

  /**
   * @param haploidSize number of different haploid hypotheses.
   * @param code encoding used for hypotheses.
   * @param hypotheses the hypotheses containing alleles to be counted.
   * @return the number of different alleles (aka haploid hypotheses) in the given set of hypotheses.
   */
  public static int numberAlleles(final int haploidSize, final Code code, final int... hypotheses) {
    //if you try and optimize this remember that the code sizes can be quite large for complex calling
    final boolean[] mark = new boolean[haploidSize];
    for (final int hyp : hypotheses) {
      mark(mark, code, hyp);
    }
    int cnt = 0;
    for (final boolean boo : mark) {
      if (boo) {
        cnt++;
      }
    }
    return cnt;
  }

  /**
   * @param haploidSize number of different haploid hypotheses.
   * @param code encoding used for hypotheses.
   * @param hyp0 a hypothesis containing alleles to be counted.
   * @param hyp1 a hypothesis containing alleles to be counted.
   * @param hyp2 a hypothesis containing alleles to be counted.
   * @return the number of different alleles (aka haploid hypotheses) in the given set of hypotheses.
   */
  public static int numberAlleles(final int haploidSize, final Code code, final int hyp0, final int hyp1, final int hyp2) {
    //if you try and optimize this remember that the code sizes can be quite large for complex calling
    final boolean[] mark = new boolean[haploidSize];
    mark(mark, code, hyp0);
    mark(mark, code, hyp1);
    mark(mark, code, hyp2);
    int cnt = 0;
    for (final boolean boo : mark) {
      if (boo) {
        cnt++;
      }
    }
    return cnt;
  }

  /**
   * @param haploidSize number of different haploid hypotheses.
   * @param code encoding used for hypotheses.
   * @param excl the hypothesis to be excluded (if -1 then exclude nothing).
   * @param hypotheses the hypotheses containing alleles to be counted.
   * @return the number of different alleles (aka haploid hypotheses) in the given set of hypotheses less a specified hypothesis.
   */
  public static int numberAllelesExclude(final int haploidSize, final Code code, final int excl, final int... hypotheses) {
    //if you try and optimize this remember that the code sizes can be quite large for complex calling
    final boolean[] mark = new boolean[haploidSize];
    for (final int hyp : hypotheses) {
      mark(mark, code, hyp);
    }
    if (excl >= 0) {
      mark[excl] = false;
    }
    int cnt = 0;
    for (final boolean boo : mark) {
      if (boo) {
        cnt++;
      }
    }
    return cnt;
  }

  /**
   * @param haploidSize number of different haploid hypotheses.
   * @param code encoding used for hypotheses.
   * @param excl the hypothesis to be excluded (if -1 then exclude nothing).
   * @param hyp0 a hypothesis containing alleles to be counted.
   * @param hyp1 a hypothesis containing alleles to be counted.
   * @param hyp2 a hypothesis containing alleles to be counted.
   * @return the number of different alleles (aka haploid hypotheses) in the given set of hypotheses less a specified hypothesis.
   */
  public static int numberAllelesExclude(final int haploidSize, final Code code, final int excl, final int hyp0, final int hyp1, final int hyp2) {
    //if you try and optimize this remember that the code sizes can be quite large for complex calling
    final boolean[] mark = new boolean[haploidSize];
    mark(mark, code, hyp0);
    mark(mark, code, hyp1);
    mark(mark, code, hyp2);
    if (excl >= 0) {
      mark[excl] = false;
    }
    int cnt = 0;
    for (final boolean boo : mark) {
      if (boo) {
        cnt++;
      }
    }
    return cnt;
  }

  /**
   * Returns whether the given family is callable, in that the output genomes have at least one child and that child has one or both parents.
   * @param outputGenomes genomes that are to be output
   * @param family a nuclear family
   * @return whether family variants can be called
   */
  public static boolean isCallableAsFamily(List<String> outputGenomes, Family family) {
    int inOutput = 0;
    for (String member : family.getMembers()) {
      if (outputGenomes.contains(member)) {
        inOutput++;
      }
    }
    int childrenInOutput = 0;
    for (String member : family.getChildren()) {
      if (outputGenomes.contains(member)) {
        childrenInOutput++;
      }
    }
    return inOutput > 1 && childrenInOutput > 0;
  }

  /**
   * Remember information from command line and reference file to enable ploidy of a sequence to be determined given the
   * sex of a sample.
   * @param params command line parameters.
   * @return the SexMemo
   * @throws java.io.IOException whenever.
   */
  public static SexMemo createSexMemo(final VariantParams params) throws IOException {
    return new SexMemo(params.genome() == null ? null :  params.genome().reader(), params.ploidy());
  }


  /**
   * Determine whether any of the models have non-reference evidence.
   * @param models list of models
   * @return true if any model has non-reference evidence
   */
  public static boolean hasOnlyRefCoverage(final List<ModelInterface<?>> models) {
    // If every model is equal to the reference no call is required.
    for (final ModelInterface<?> m : models) {
      if (m.statistics().nonRefCount() > 0) {
        return false;
      }
    }
    return true;
  }

  /**
   * Calculate the mean ratio of ambiguous reads
   * @param models list of fully incremented models
   * @return the mean of the ambiguities of each model
   */
  public static double meanAmbiguityRatio(final List<ModelInterface<?>> models) {
    double ambig = 0;
    int acount = 0;
    for (final ModelInterface<?> b : models) {
      final Double aa = b.statistics().ambiguityRatio();
      if (aa != null) {
        ambig += aa;
        acount++;
      }
    }
    //    System.err.println("" + " ambig/models.size()=" + (ambig/models.size()) );
    return acount == 0 ? 0 : ambig / acount;
  }

  /**
   * Gets the total coverage of any of the supplied models
   * @param models list of models
   * @return the total coverage
   */
  public static int totalCoverage(final List<ModelInterface<?>> models) {
    int coverage = 0;
    for (final ModelInterface<?> b : models) {
      coverage += b.statistics().coverage();
    }
    return coverage;
  }

  /**
   * Gets the maximum coverage of any of the supplied models
   * @param models list of models
   * @return the highest coverage
   */
  public static int maxCoverage(final List<ModelInterface<?>> models) {
    int coverage = 0;
    for (final ModelInterface<?> b : models) {
      coverage = Math.max(coverage, b.statistics().coverage());
    }
    return coverage;
  }
}
