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
package com.rtg.index.hash.ngs.instances;


import static com.rtg.util.StringUtils.LS;

import java.io.IOException;

import com.rtg.index.hash.ngs.HashFunctionFactory;
import com.rtg.index.hash.ngs.NgsHashFunction;
import com.rtg.index.hash.ngs.instances.AbstractSplitTest.ReadCallAccumulate;
import com.rtg.index.hash.ngs.instances.AbstractSplitTest.TemplateCallCheck;

import org.junit.Assert;

/**
 * Check that all possible substitutions in a string are found when filtered via the <code>HashFunction</code>.
 */
public class SubstituteCG {
  private static final int PADDING = 8;
  private final String mString;
  private final int mLength;
  private final HashFunctionFactory mHF;
  private final boolean mCheckFound;
  private final int mGap;
  private final int mOverlap;

  private int mNotFound;
  private int mTotal;
  private StringBuilder mResults = new StringBuilder("not run");
  private final char[] mChars;
  private final char[] mSubstitutions;

  /**
   * @param string to be used for doing substitution test.
   * @param factory factory for generating mask.
   * @param checkFound if true check that everything is actually found
   * @param cgGap size of right gap (version 1).
   * @param cgOverlap size of overlap (version 1 and version 2)
   */
  public SubstituteCG(final String string, final HashFunctionFactory factory, final boolean checkFound, final int cgGap, int cgOverlap) {
    mString = string;
    mLength = mString.length();
    mChars = new char[mLength];
    mHF = factory;
    mCheckFound = checkFound;
    mGap = cgGap;
    mOverlap = cgOverlap;
    mSubstitutions = new char[mLength];
  }

  /**
   * @param error Maximum number of errors which will be generated.
   * @return the number of errors that were not found
   * @throws IOException If an I/O error occurs
   */
  public int substituteProtected(final int error) throws IOException {
    mNotFound = 0;
    mTotal = 0;
    mResults = new StringBuilder();
    substitutePrivate(error, 0);
    return mNotFound;
  }

  @Override
  public String toString() {
    return "Missed " + mNotFound + "/" + mTotal;
  }

  /**
   * @return textual summary of how many combinations were not found
   */
  public String details() {
    return mResults.toString();
  }

  static String repeat(final String rep, final int n) {
    final StringBuilder sb = new StringBuilder();
    for (int i = 0; i < n; ++i) {
      sb.append(rep);
    }
    return sb.toString();
  }

  static String cgToTemplate(final String mutant, final int gap, int overlap) {
    if (29 == mutant.length()) {
      return repeat("a", PADDING)
        + mutant.substring(0, 10)
        + mutant.substring(10 + overlap, mutant.length())
        + repeat("a", PADDING);
    } else if (35 == mutant.length()) {
      final StringBuilder sb = new StringBuilder(repeat("a", PADDING));
      sb.append(mutant.substring(0, 5));
      int pos = 5 + overlap;
      sb.append(mutant.substring(pos, pos + 20 - overlap)); // Dodgy
      pos += 20 - overlap;
      sb.append(repeat("a", gap));
      sb.append(mutant.substring(pos));
      sb.append(repeat("a", PADDING));
      return sb.toString();
    } else {
      Assert.fail();
      return null;
    }
  }

  private void substitutePrivate(final int error, final int soFar) throws IOException {
    //System.err.println("error=" + error + " soFar=" + soFar);
    if (error == 0 || soFar == mLength) {
      for (int i = soFar; i < mLength; ++i) {
        mChars[i] = mString.charAt(i);
      }
      final ReadCallAccumulate rc = new ReadCallAccumulate();
      final TemplateCallCheck tcc = new TemplateCallCheck(rc.map());
      final NgsHashFunction hf = mHF.create(rc, tcc);
      AbstractSplitTest.encode(hf, mString);
      hf.readAll(0, false);
      hf.reset();
      final String mutant = new String(mChars);
      final String muttie = cgToTemplate(mutant, mGap, mOverlap);
      for (int i = 0; i < muttie.length(); ++i) {
        AbstractSplitTest.encode(hf, muttie, i);
        hf.templateForward(0);
      }
      if (mCheckFound && !tcc.mFoundDone) {
        mResults.append("R:" + mString + LS);
        mResults.append("T:" + muttie.substring(PADDING, muttie.length() - PADDING) + LS);
        mResults.append("  " + new String(mSubstitutions, 0, soFar) + LS);
        ++mNotFound;
      }
      ++mTotal;
      return;
    }
    //make one substitution and as well the non-substitution case
    for (int j = 1, c = 0; j < AbstractSplitTest.CHARS.length; ++j) {
      mChars[soFar]  = AbstractSplitTest.CHARS[j];
      if (mChars[soFar] == mString.charAt(soFar)) {
        mSubstitutions[soFar] = ("" + (soFar % 10)).charAt(0);
        substitutePrivate(error, soFar + 1);
      } else {
        if (c == 0 && error > 0) {
          mSubstitutions[soFar] = '.';
          substitutePrivate(error - 1, soFar + 1);
        }
        ++c;
      }
    }
  }
}

