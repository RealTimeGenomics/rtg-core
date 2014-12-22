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

import java.util.ArrayList;
import java.util.List;

import com.rtg.util.License;

/**
 * Enable easy reporting of such things as counters to the developer log.
 */
public final class Spy {

  private static final List<Object> SPIES = new ArrayList<>();

  /**
   * Add a new spy.
   * @param spy to be added.
   */
  public static void add(Object spy) {
    SPIES.add(spy);
  }

  /**
   * Generate a report to the log from all the current spies.
   */
  public static void report() {
    //System.err.println("Spy report:" + SPIES.size());
    if (License.isDeveloper()) {
      for (final Object obj : SPIES) {
        Diagnostic.developerLog(obj.toString());
      }
    }
  }

  private Spy() { }
}
