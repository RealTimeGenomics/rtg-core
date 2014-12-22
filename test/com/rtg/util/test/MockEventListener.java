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

package com.rtg.util.test;

import com.rtg.util.diagnostic.DiagnosticEvent;
import com.rtg.util.diagnostic.DiagnosticListener;
import com.rtg.util.diagnostic.ErrorEvent;

import junit.framework.Assert;

/**
 * Make it easy to test Diagnostic messages.
 *
 */
public final class MockEventListener implements DiagnosticListener {

  DiagnosticEvent<?> mEvent = null;

  @Override
  public void handleDiagnosticEvent(final DiagnosticEvent<?> event) {
    mEvent = event;
  }

  /**
   * Check that an error event was recieved with the supplied message.
   * @param expected the expected message
   * @return true if the message was recieved and contained the message.
   */
  public boolean compareErrorMessage(final String expected) {
    Assert.assertNotNull("Event is null", mEvent);
    boolean res;
    res = mEvent instanceof ErrorEvent;
    //System.err.println(mEvent.getMessage());
    res &= expected.equals(mEvent.getMessage());
    return res;
  }

  /**
   * @return the event recieved.
   */
  public DiagnosticEvent<?> getEvent() {
    return mEvent;
  }

  /**
   */
  @Override
  public void close() {
  }
}
