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

import java.io.EOFException;
import java.io.IOException;

import junit.framework.TestCase;

/**
 * Test class
 */
public class IORunnableProxyTest extends TestCase {

  public void testSomeMethod() throws Exception {
    final IORunnableProxy proxIO = new IORunnableProxy(new IORunnable() {
      @Override
      public void run() throws IOException {
        throw new EOFException("generic IO exception");
      }
    });
    Thread t = new Thread(proxIO);
    t.start();
    t.join();
    try {
      proxIO.checkError();
      fail();
    } catch (EOFException e) {
    }

    final IORunnableProxy proxRuntime = new IORunnableProxy(new IORunnable() {
      @Override
      public void run() {
        throw new IllegalArgumentException("gotcha");
      }
    });
    t = new Thread(proxRuntime);
    t.start();
    t.join();
    try {
      proxRuntime.checkError();
      fail();
    } catch (IllegalArgumentException e) {
    }

    final IORunnableProxy proxError = new IORunnableProxy(new IORunnable() {
      @Override
      public void run() {
        throw new OutOfMemoryError("not really");
      }
    });
    t = new Thread(proxError);
    t.start();
    t.join();
    try {
      proxError.checkError();
      fail();
    } catch (OutOfMemoryError e) {
    }

    final IORunnableProxy proxOk = new IORunnableProxy(new IORunnable() {
      @Override
      public void run() {
      }
    });
    t = new Thread(proxOk);
    t.start();
    t.join();
    proxOk.checkError();
  }

}
