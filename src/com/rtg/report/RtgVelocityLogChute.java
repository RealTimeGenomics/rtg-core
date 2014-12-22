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
package com.rtg.report;

import org.apache.velocity.runtime.RuntimeServices;
import org.apache.velocity.runtime.log.LogChute;

import com.rtg.util.License;
import com.rtg.util.diagnostic.Diagnostic;

/**
 * For Stupid velocity
 */
public class RtgVelocityLogChute implements LogChute {
  @Override
  public void init(RuntimeServices rs) {
  }

  @Override
  public void log(int level, String message) {
    if (isLevelEnabled(level)) {
      if (level >= INFO_ID) {
        Diagnostic.userLog(message);
      } else if (level == DEBUG_ID) {
        Diagnostic.developerLog(message);
      }
    }
  }

  @Override
  public void log(int level, String message, Throwable t) {
    if (isLevelEnabled(level)) {
      if (level >= INFO_ID) {
        Diagnostic.userLog(message);
      } else if (level == DEBUG_ID) {
        Diagnostic.developerLog(message);
      }
    }
  }

  @Override
  public boolean isLevelEnabled(int level) {
    return License.isDeveloper() && level >= DEBUG_ID || level >= INFO_ID;
  }
}
