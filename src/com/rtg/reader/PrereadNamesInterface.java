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

package com.rtg.reader;

import java.io.IOException;
import java.io.OutputStream;

/**
 * Container for sequence names accessed by sequence id
 */
public interface PrereadNamesInterface {

  /**
   * Returns number of names
   * @return length of names
   */
  long length();

  /**
   * Returns the sequence name for a specified id.
   * @param id sequence id
   * @return sequence name
   */
  String name(final long id);

  /**
   * Calculate the checksum of the names in a manner compatible with
   * how the checksum is calculated in the SDF.
   *
   * @return the checksum of the names.
   */
  long calcChecksum();

  /**
   * @return size of object in bytes
   */
  long bytes();

  /**
   * Convenience method to append a name to an Appendable.  This avoid
   * the string creation of the <code>name()</code> method.
   *
   * @param a an <code>Appendable</code> value
   * @param id a <code>long</code> value
   * @throws IOException when writing to the appendable.
   */
  void writeName(final Appendable a, final long id) throws IOException;

  /**
   * Convenience method to append a name to an Appendable.  This avoid
   * the string creation of the <code>name()</code> method.
   *
   * @param os an <code>OutputStream</code> value
   * @param id a <code>long</code> value
   * @throws IOException when writing to the stream.
   */
  void writeName(OutputStream os, long id) throws IOException;

}
