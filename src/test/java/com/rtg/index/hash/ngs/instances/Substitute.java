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


import java.io.IOException;

import com.rtg.index.hash.ngs.HashFunctionFactory;
import com.rtg.index.hash.ngs.NgsHashFunction;
import com.rtg.index.hash.ngs.instances.AbstractSplitTest.ReadCallAccumulate;
import com.rtg.index.hash.ngs.instances.AbstractSplitTest.TemplateCallCheck;

import org.junit.Assert;

/**
 * Check that all possible substitutions in a string are found when filtered via the <code>HashFunction</code>.
 */
public class Substitute {
  private final String mString;
  private final int mLength;
  private final HashFunctionFactory mHF;
  private final boolean mCheckFound;

  private boolean mNoErrors;
  private final char[] mChars;
  private final char[] mSubstitutions;

  /**
   * @param string to be used for doing substitution test.
   * @param factory factory for generating mask.
   * @param checkFound if true check that everything is actually found
   */
  public Substitute(final String string, final HashFunctionFactory factory, final boolean checkFound) {
    mString = string;
    mLength = mString.length();
    mChars = new char[mLength];
    mHF = factory;
    mCheckFound = checkFound;
    mSubstitutions = new char[mLength];
  }

  /**
   *
   * @param error Maximum number of errors which will be generated.
   * @throws IOException If an I/O error occurs
   */
  public void substituteProtected(final int error) throws IOException {
    mNoErrors = true;
    substitutePrivate(error, 0);
    Assert.assertTrue(mNoErrors);
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
      AbstractSplitTest.encode(hf, mutant);
      hf.templateForward(0);
      if (mCheckFound && !tcc.mFoundDone) {
        System.err.println("errorFound mLength: " + mLength + " soFar: " + soFar + " error:" + error);
        System.err.println(new String(mSubstitutions, 0, soFar));
        mNoErrors = false;
      }
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

