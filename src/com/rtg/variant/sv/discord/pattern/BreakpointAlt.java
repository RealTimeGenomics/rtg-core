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

package com.rtg.variant.sv.discord.pattern;

/**
 * Class that works out what kind of breakpoint a VCF alt record represents
*         Date: 15/03/12
*         Time: 3:29 PM
*/
class BreakpointAlt {
  final boolean mLocalUp;
  final boolean mRemoteUp;
  final String mRemoteChr;
  final int mRemotePos; // 1-based
  private static final String UP_BRACKET = "]";
  private static final String DOWN_BRACKET = "[";

  BreakpointAlt(String alt) {
    final int tmpPos;
    final int bracketPos;
    final String bracket;
//    System.err.println(alt.indexOf(UP_BRACKET));
    if ((tmpPos = alt.indexOf(UP_BRACKET)) > -1) {
      mRemoteUp = true;
      bracket = UP_BRACKET;
      bracketPos = tmpPos;
    } else {
      bracket = DOWN_BRACKET;
      bracketPos = alt.indexOf(DOWN_BRACKET);
      mRemoteUp = false;
    }
//    System.err.println("" + " bracketPos=" + bracketPos + " bracket=" + bracket );
    mLocalUp  = !(bracketPos == 0);
    final String remoteString = alt.substring(bracketPos + 1, alt.indexOf(bracket, bracketPos + 1));
    mRemoteChr = remoteString.substring(0, remoteString.indexOf(":"));
    mRemotePos = Integer.parseInt(remoteString.substring(remoteString.indexOf(":") + 1));
  }

  /**
   * @return True if the remote part of the breakpoint is orientation up see {@link com.rtg.variant.sv.discord.Orientation}
   */
  public boolean isRemoteUp() {
    return mRemoteUp;
  }

  /**
   * @return the name of the chromosome this breakpoint joins to
   */
  public String getRemoteChr() {
    return mRemoteChr;
  }

  /**
   * @return the position within the remote chromosome this breakpoint is at, one-based
   */
  public int getRemotePos() {
    return mRemotePos;
  }

  /**
   * @return True if the local part of the breakpoint is orientation up see {@link com.rtg.variant.sv.discord.Orientation}
   */
  public boolean isLocalUp() {
    return mLocalUp;
  }

}
