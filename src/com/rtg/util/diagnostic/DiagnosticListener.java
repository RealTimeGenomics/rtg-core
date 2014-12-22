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
package com.rtg.util.diagnostic;

/**
 * Interface for classes wanting to receive SLIM diagnostics.
 *
 */
public interface DiagnosticListener {

  /**
   * Called when a diagnostic is produced.  The event can be ignored if desired.
   * Implementations should not modify the passed event in any way.  The null
   * event will never be passed as input to this function.
   *
   * @param event the event
   */
  void handleDiagnosticEvent(DiagnosticEvent<?> event);

  /**
   * Allow the listener to perform any desired clean up.
   */
  void close();
}

