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
package com.rtg.util;

import static com.rtg.util.StringUtils.FS;
import static com.rtg.util.StringUtils.LS;

import java.io.IOException;
import java.io.OutputStream;

/**
 * byte and byte[] utilities.
 * Intended to help people writing directly to <code>OutputStream</code> s
 *
 */
public final class ByteUtils {

  /** System dependent line separator. */
  private static final byte[] LS_BYTES = LS.getBytes();
  /**  System dependent file separator. */
  public static final byte FS_BYTE = FS.getBytes()[0];
  /** Tab. */
  public static final byte TAB_BYTE = (byte) '\t';
  /** System independent newline separator */
  private static final byte[] NEWLINE_BYTES = "\n".getBytes();

  /**
   * Write an end of line to out.
   * @param out where to put the output.
   * @throws IOException from out.
   */
  public static void writeLn(final OutputStream out) throws IOException {
    out.write(LS_BYTES);
  }

  /**
   * Write a newline (<code>/n</code>) to out
   * @param out where to write to
   * @throws IOException from out
   */
  public static void writeNewline(final OutputStream out) throws IOException {
    out.write(NEWLINE_BYTES);
  }

  private ByteUtils() { }

}
