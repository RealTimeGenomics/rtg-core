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

import java.io.ByteArrayInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.security.DigestInputStream;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;

/**
 */
public final class MD5Utils {

  private MD5Utils() {
  }

  /**
   * Read all of the input stream and return it's MD5 digest.
   * @param in input stream to be read
   * @return MD5 digest.
   * @throws IOException when reading input stream.
   */
  public static String md5(final InputStream in) throws IOException {
    return md5(in, 1024);
  }

  /**
   * Read all of the input stream and return it's MD5 digest.
   * @param in input stream to be read
   * @return MD5 digest.
   * @throws IOException when reading input stream.
   */
  static String md5(final InputStream in, final int bufferLength) throws IOException {
    final byte[] buffer = new byte[bufferLength];
    try {
      final MessageDigest md = MessageDigest.getInstance("MD5");
      final DigestInputStream dis = new DigestInputStream(in, md);
      while (dis.read(buffer) != -1) { }
      final byte[] digest = md.digest();
      return digestToString(digest);
    } catch (final NoSuchAlgorithmException e) {
      throw new RuntimeException(e);
    }
  }

  /**
   * @param digest MD5 digest as a byte array.
   * @return MD5 digest as a string.
   */
  public static String digestToString(final byte[] digest) {
    final StringBuilder sb = new StringBuilder();
    for (byte aDigest : digest) {
      sb.append(Integer.toHexString(0x100 + (aDigest & 0xFF)).substring(1));
    }
    return sb.toString();
  }

  /**
   * Returns the MD5 for the given <code>text</code>.
   *
   * @param text a <code>String</code>
   * @return MD5
   */
  public static String md5(final String text) {
    if (text == null) {
      throw new IllegalArgumentException("String was null");
    } else if (text.length() == 0) {
      throw new IllegalArgumentException("String was 0 length");
    }
    try {
      return md5(new ByteArrayInputStream(text.getBytes()));
    } catch (IOException ioe) {
      throw new RuntimeException(ioe);  //should never be able to happen, though
    }
  }
}
