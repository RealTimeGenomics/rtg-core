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
 */
class BreakpointAlt {

  private static final char UP_BRACKET = ']';
  private static final char DOWN_BRACKET = '[';

  final boolean mLocalUp;
  final boolean mRemoteUp;
  final String mRemoteChr;
  final int mRemotePos;
  final String mRefSubs;

  BreakpointAlt(String alt) {
    final int bracketPos;
    final char bracket;
    int tmpPos;
    if ((tmpPos = alt.indexOf(UP_BRACKET)) > -1) {
      bracket = UP_BRACKET;
      bracketPos = tmpPos;
      mRemoteUp = true;
    } else if ((tmpPos = alt.indexOf(DOWN_BRACKET)) > -1) {
      bracket = DOWN_BRACKET;
      bracketPos = tmpPos;
      mRemoteUp = false;
    } else {
      throw new IllegalArgumentException("Invalid break end specification: " + alt);
    }
    mLocalUp = bracketPos != 0;
    final int nextBracketPos = alt.indexOf(bracket, bracketPos + 1);
    if (nextBracketPos == -1) {
      throw new IllegalArgumentException("Invalid break end specification: " + alt);
    }
    final String remoteString = alt.substring(bracketPos + 1, nextBracketPos);
    final int colonPos = remoteString.indexOf(':');
    if (colonPos == -1) {
      throw new IllegalArgumentException("Invalid break end specification: " + alt);
    }
    mRemoteChr = remoteString.substring(0, colonPos);
    if (mRemoteChr.indexOf('<') != -1) {
      throw new IllegalArgumentException("Contig break ends are not supported: " + alt); // e.g. "C[<ctg1>:7["
    }
    mRemotePos = Integer.parseInt(remoteString.substring(colonPos + 1)) - 1; // to 0-based
    mRefSubs = mLocalUp ? alt.substring(0, bracketPos) : alt.substring(nextBracketPos + 1);
  }

  BreakpointAlt(String refSubs, boolean localUp, String remoteChr, int remotePos, boolean remoteUp) {
    mRefSubs = refSubs;
    mLocalUp = localUp;
    mRemoteChr = remoteChr;
    mRemotePos = remotePos;
    mRemoteUp = remoteUp;
  }

  /**
   * @return true if the remote part of the breakpoint is orientation up see {@link com.rtg.variant.sv.discord.Orientation}
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
   * @return true if the local part of the breakpoint is orientation up see {@link com.rtg.variant.sv.discord.Orientation}
   */
  public boolean isLocalUp() {
    return mLocalUp;
  }

  /**
   * @return the bases that substitute for the reference bases
   */
  public String getRefSubstitution() {
    return mRefSubs;
  }

  @Override
  public String toString() {
    final char bracket = mRemoteUp ? UP_BRACKET : DOWN_BRACKET;
    return mLocalUp
      ? mRefSubs + bracket + mRemoteChr + ":" + (mRemotePos + 1) + bracket
      : bracket + mRemoteChr + ":" + (mRemotePos + 1) + bracket + mRefSubs ;
  }

}
