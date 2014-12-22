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

package com.rtg.util.store;

import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;

/**
 */
public interface StoreFile {

  /**
   * @return an output stream to write to this file.
   * @throws IOException if problem creating the output stream.
   */
  OutputStream outputStream() throws IOException;

  /**
   * @return an input stream from this file.
   * @throws IOException  if problem creating the input stream.
   */
  InputStream inputStream() throws IOException;

  /**
   * @return the contents of this file.
   * @throws IOException if problems reading the contents.
   */
  String content() throws IOException;

  /**
   * @return the name of this child.
   */
  String name();
}
