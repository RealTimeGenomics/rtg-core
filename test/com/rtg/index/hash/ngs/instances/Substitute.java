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
      for (int i = soFar; i < mLength; i++) {
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
    for (int j = 1, c = 0; j < AbstractSplitTest.CHARS.length; j++) {
      mChars[soFar]  = AbstractSplitTest.CHARS[j];
      if (mChars[soFar] == mString.charAt(soFar)) {
        mSubstitutions[soFar] = ("" + (soFar % 10)).charAt(0);
        substitutePrivate(error, soFar + 1);
      } else {
        if (c == 0 && error > 0) {
          mSubstitutions[soFar] = '.';
          substitutePrivate(error - 1, soFar + 1);
        }
        c++;
      }
    }
  }
}

