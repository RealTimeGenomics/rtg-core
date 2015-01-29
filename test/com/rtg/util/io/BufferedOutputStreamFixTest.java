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
package com.rtg.util.io;

import java.io.BufferedOutputStream;
import java.io.IOException;
import java.io.OutputStream;

import junit.framework.TestCase;

/**
 */
public class BufferedOutputStreamFixTest extends TestCase {

  private static final String FLUSH_MSG = "I don't like to flush";

  private class BadOutputStream extends OutputStream {

    boolean mBadstate = true;

    @Override
    public void write(int b) {
    }

    @Override
    public void flush() throws IOException {
      if (mBadstate) {
        mBadstate = false;
        throw new IOException(FLUSH_MSG);
      }
    }

    @Override
    public void close() throws IOException {
      throw new IOException("I don't like to close");
    }
  }

  public void test() {
    //Java7 contains this bug, and it appears to have been fixed in Java8.
    //first test if the bug still exists in the java being run
//    try {
//      try (final BufferedOutputStream bos = new BufferedOutputStream(new BadOutputStream())) {
//        bos.write("token".getBytes());
//      }
//    } catch (IOException e) {
//      if (FLUSH_MSG.equals(e.getMessage())) {
//        System.err.println("Seems java has now fixed their bug in FilterOutputStream");
//      }
//    }

    //then check we fixed it
    try {
      try (final BufferedOutputStream bos = new BufferedOutputStreamFix(new BadOutputStream())) {
        bos.write("token".getBytes());
      }
      fail("Should have thrown exception");
    } catch (IOException e) {
      assertEquals(FLUSH_MSG, e.getMessage());
    }
  }
}
