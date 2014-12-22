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
package com.rtg.util;

import java.io.IOException;


/**
 * Throw IOException on all operations.
 */
public class MockAppendable implements Appendable {
  /**
   */
  @Override
  public Appendable append(final char arg0) throws IOException {
    throw new IOException();
  }
  /**
   */
  @Override
  public Appendable append(final CharSequence arg0, final int arg1, final int arg2)
  throws IOException {
    throw new IOException();
  }
  /**
   */
  @Override
  public Appendable append(final CharSequence arg0) throws IOException {
    throw new IOException();
  }

}
