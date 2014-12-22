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
package com.rtg.index.hash.ngs;

import java.io.IOException;

/**
 * Common functions to both <code>ReadHashFunction</code> and <code>TemplateHashFunction</code>.
 */
public interface ReadHashFunction extends CommonHashFunction {


  /**
   * Force creation of an array to hold the read bit sequences from the <code>HashLoop</code> supplied to the <code>HashFunction</code>
   * @param numberReads number of read sequences.
   */
  void setReadSequences(long numberReads);

  /**
   * Process all windows for the specified read.
   * @param readId number of the read (&gt;=0).
   * @param reverse if true store reverse complement.
   * @throws IOException If an I/O error occurs
   */
  void readAll(final int readId, final boolean reverse) throws IOException;

  /**
   * Store the current values in <code>readSequences</code>.
   * @param id2 read id.
   * @param reverse if true use reverse complement.
   */
  void setValues(int id2, boolean reverse);
}

