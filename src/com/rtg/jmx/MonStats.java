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
package com.rtg.jmx;

import java.io.IOException;

/**
 * Interface for objects adding monitor stats to output.
 */
public interface MonStats {

  /**
   * Output any data to be produced once at the start of monitoring.
   *
   * @param out the destination <code>Appendable</code>
   * @throws IOException if there is a problem during output
   */
  void addHeader(Appendable out) throws IOException;

  /**
   * Output tab separated column labels corresponding to data produced by this object.
   *
   * @param out the destination <code>Appendable</code>
   * @throws IOException if there is a problem during output
   */
  void addColumnLabelsTop(Appendable out) throws IOException;

  /**
   * Output tab separated column labels corresponding to data produced by this object.
   *
   * @param out the destination <code>Appendable</code>
   * @throws IOException if there is a problem during output
   */
  void addColumnLabelsBottom(Appendable out) throws IOException;

  /**
   * Output tab separated column data.
   *
   * @param out the destination <code>Appendable</code>
   * @throws IOException if there is a problem during output
   */
  void addColumnData(Appendable out) throws IOException;

}
