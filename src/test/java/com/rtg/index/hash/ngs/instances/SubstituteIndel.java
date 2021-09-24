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
package com.rtg.index.hash.ngs.instances;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;

import com.rtg.index.hash.ngs.HashFunctionFactory;
import com.rtg.index.hash.ngs.NgsHashFunction;
import com.rtg.index.hash.ngs.instances.AbstractSplitTest.ReadCallAccumulate;
import com.rtg.index.hash.ngs.instances.AbstractSplitTest.TemplateCallCheck;

/**
 * Check that all possible substitutions in a string are found when filtered via the <code>HashFunction</code>.
 */
public class SubstituteIndel {

  /**
   * Check all mutations including indels.
   * @param factory for creating mask.
   * @param str to be mutated.
   * @param error number of indels allowed.
   * @param cgGap number of nt to insert in gap for CG data.
   * @throws IOException when ever.
   */
  public static void checkIndel(final HashFunctionFactory factory, final String str, final int error, final int cgGap) throws IOException {
    final int length = str.length();
    final SubstituteIndel sub = new SubstituteIndel(str, error);
    final Collection<String> templates = sub.substitute(error);
    final String pref = "cta".substring(0, error);
    final String suff = "gct".substring(0, error);
    final ReadCallAccumulate rc = new ReadCallAccumulate();
    final TemplateCallCheck tcc = new TemplateCallCheck(rc.map());
    final NgsHashFunction hf = factory.create(rc, tcc);
    for (final String templ : templates) {
      //put padding on the end to allow up to 3 indels
      AbstractSplitTest.encode(hf, str);
      hf.readAll(0, false);
      hf.reset();
      final String template = pref + templ + suff;
      for (int i = 0; i + length <= template.length(); ++i) {
        final String mutant = template.substring(i, i + length);
        final String muttie = cgGap > 0 ? SubstituteCG.cgToTemplate(mutant, cgGap, 0) : mutant;
        for (int j = 0; j < muttie.length(); ++j) {
          AbstractSplitTest.encode(hf, muttie, j);
          hf.templateForward(0);
        }
        hf.reset();
      }
      if (!tcc.mFound) {
        System.err.println(str);
        System.err.println(template);
        AbstractSplitTest.fail();
      }
    }
  }

  private final String mString;
  private final int mLength;
  private char[] mChars = null;
  private final ArrayList<String> mResult = new ArrayList<>();

  SubstituteIndel(final String string, final int error) {
    mString = string;
    mLength = mString.length();
    assert error >= 0 && error <= 3;
  }

  Collection<String> substitute(final int error) {
    assert error >= 0 && error <= 3;
    mChars = new char[mLength + 2 * error];
    substitute(error, 0, 0);
    //    System.err.println("base count=" + mStats0);
    //    System.err.println("total count=" + mStats1);
    return mResult;
  }

  private void substitute(final int error, final int soFar, final int out) {
    //System.err.println("error=" + error + " soFar=" + soFar);
    if (error == 0 || soFar == mLength || out == mLength) {
      int j = out;
      for (int i = soFar; j < mLength || i < mLength; ++i, ++j) {
        if (i >= mLength) {
          mChars[j] = 't';
        } else {
          mChars[j] = mString.charAt(i);
        }
      }
      mResult.add(new String(mChars).substring(0, j));
      return;
    }
    //delete
    substitute(error - 1, soFar + 1, out);
    //insert
    mChars[out]  = AbstractSplitTest.CHARS[2];
    substitute(error - 1, soFar, out + 1);
    //make one substitution and as well the non-substitution case
    for (int j = 1, c = 0; j < AbstractSplitTest.CHARS.length; ++j) {
      mChars[out]  = AbstractSplitTest.CHARS[j];
      if (mChars[out] == mString.charAt(soFar)) {
        substitute(error, soFar + 1, out + 1);
      } else {
        if (c == 0 && error > 0) {
          substitute(error - 1, soFar + 1, out + 1);
        }
        ++c;
      }
    }
  }
}

