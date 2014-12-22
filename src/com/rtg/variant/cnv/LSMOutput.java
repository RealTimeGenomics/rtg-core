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
package com.rtg.variant.cnv;

import java.io.IOException;

/**
 */
public interface LSMOutput {

  /**
   * Perform an output action for each fitted region.
   * @param id current sequence id.
   * @param start start position within the sequence.
   * @param end end position within the sequence.
   * @param sum sigma_i=start^end x_i
   * @param sum2 sigma_i=start^end x_i^2
   * @throws IOException IO Exceptions happen sometimes
   */
  void out(final String id, final int start, final int end, final double sum, final double sum2) throws IOException;

  /**
   * Close any files or other resources.
   * @throws IOException IO Exceptions happen sometimes
   */
  void close() throws IOException;
}
