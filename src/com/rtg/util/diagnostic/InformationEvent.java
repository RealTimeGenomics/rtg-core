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
 * Event used to report informational messages.
 *
 */
public final class InformationEvent extends DiagnosticEvent<InformationType> {

  /**
   * Constructs an information event.
   * @param type the type of the event
   * @param params the parameters for the event (type specific)
   */
  public InformationEvent(final InformationType type, final String... params) {
    super(type, params);
  }

}
