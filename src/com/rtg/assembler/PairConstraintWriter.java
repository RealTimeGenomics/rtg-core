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

package com.rtg.assembler;

import java.io.Closeable;
import java.io.PrintStream;

/**
 */
public class PairConstraintWriter implements Closeable {
  final PrintStream mOut;
  final ConstraintCache mCache = new ConstraintCache();


  PairConstraintWriter(PrintStream out) {
    mOut = out;
  }

  void writeConstraint(long contig1, long contig2, int distanceFromEnd1, int distanceFromEnd2, boolean leftEnd, boolean rightEnd) {
    mOut.println(String.valueOf(contig1) + "\t" + contig2 + "\t" + distanceFromEnd1 + "\t" + distanceFromEnd2 + "\t" + leftEnd + "\t" + rightEnd);
  }
  @Override
  public void close() {
    mOut.close();
  }

}
