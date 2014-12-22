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
 * Wrapper for running an {@link IORunnable} inside a framework that only accepts
 * {@link Runnable}
 */
public class IORunnableProxy implements Runnable {

  private final IORunnable mWrapped;

  private IOException mIOException;
  private RuntimeException mRuntimeException;
  private Error mError;

  /**
   * @param wrapped task to wrap
   */
  public IORunnableProxy(IORunnable wrapped) {
    mWrapped = wrapped;
  }


  @Override
  public void run() {
    try {
      mWrapped.run();
    } catch (IOException e) {
      mIOException = e;
    } catch (RuntimeException e) {
      mRuntimeException = e;
    } catch (Error e) {
      mError = e;
    }
  }

  /**
   * Throws any exception or error encountered during the run, wrapped in the appropriate wrap. If there were
   * no problems nothing happens
   * @throws IOException if an <code>IOException</code> was encountered during run
   * @throws RuntimeException if a <code>RuntimeException</code> was encountered during run.
   * @throws Error if a <code>Error</code> was encountered during run.
   */
  public void checkError() throws IOException {
    if (mIOException != null) {
      throw mIOException;
    } else if (mRuntimeException != null) {
      throw mRuntimeException;
    } else if (mError != null) {
      throw mError;
    }
  }
}
