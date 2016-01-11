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

package com.rtg.metagenomics.metasnp;

import java.io.IOException;
import java.util.List;

/**
 * Source of variants for metasnp.
 */
public interface MetaSnpReader extends AutoCloseable {

  /**
   * @return the list of sample names from the header or null if the header is missing
   */
  List<String> samples();

  /**
   * @return next line of input
   * @throws IOException if an I/O error occurs
   */
  MetaSnpLine nextLine() throws IOException;

  @Override
  void close() throws IOException;
}
