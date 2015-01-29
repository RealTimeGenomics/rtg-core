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
package com.rtg.util.io;


import java.io.ByteArrayOutputStream;
import java.io.EOFException;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.Reader;
import java.io.StringWriter;
import java.net.URL;

/**
 * Access to io functionality.
 *
 */
public final class IOUtils {

  private IOUtils() { }

  /** Length of buffer to use during reading */
  private static final int BUFFER_LENGTH = 16384;

  private static final int EOF = -1;

  /**
   * Reads all the data from the supplied URL into a byte array.
   *
   * @param url the URL
   * @return a byte array containing the stream data.
   * @exception IOException if an error occurs during IO.
   */
  public static byte[] readData(final URL url) throws IOException {
    try (InputStream input = url.openStream()) {
      return readData(input);
    }
  }

  /**
   * Reads all the data from the supplied file into a byte array.
   *
   * @param input the file
   * @return a byte array containing the stream data.
   * @exception IOException if an error occurs during IO.
   */
  public static byte[] readData(File input) throws IOException {
    try (InputStream is = new FileInputStream(input)) {
      return readData(is);
    }
  }

  /**
   * Reads all the data from the supplied InputStream into a byte array.
   *
   * @param input the InputStream
   * @return a byte array containing the stream data.
   * @exception IOException if an error occurs during IO.
   */
  public static byte[] readData(final InputStream input) throws IOException {
    final byte[] inputBuffer = new byte[BUFFER_LENGTH];
    final ByteArrayOutputStream byteOutput = new ByteArrayOutputStream(BUFFER_LENGTH);
    int bytesRead;
    while ((bytesRead = input.read(inputBuffer)) != -1) {
      byteOutput.write(inputBuffer, 0, bytesRead);
    }
    final byte[] r = byteOutput.toByteArray();
    byteOutput.close();
    return r;
  }

  /**
   * Reads specified amount of data from stream or amount left in stream.
   * @param is stream to read from
   * @param buf buffer to read into
   * @param offset start position in buffer to read into
   * @param length amount of data to read in bytes
   * @return amount of data read
   * @throws IOException if an IO error occurs
   */
  public static int readAmount(InputStream is, byte[] buf, int offset, int length) throws IOException {
    int read = 0;
    int len;
    while (read < length && (len = is.read(buf, offset + read, length - read)) > 0) {
      read += len;
    }
    return read;
  }

  /**
   * Reads specified amount of data from stream.
   * @param is stream to read from
   * @param buf buffer to read into
   * @param offset start position in buffer to read into
   * @param length amount of data to read in bytes
   * @throws IOException If an IO error occurs
   */
  public static void readFully(InputStream is, byte[] buf, int offset, int length) throws IOException {
    final int read = readAmount(is, buf, offset, length);
    if (read < length) {
      throw new EOFException();
    }
  }

  /**
   * Read all of a URL into a String.
   *
   * @param url the URL to read.
   * @return a String containing the contents of the URL
   * @exception IOException If there is a problem during reading.
   */
  public static String readAll(final URL url) throws IOException {
    try (InputStream input = url.openStream()) {
      return readAll(input);
    }
  }


  /**
   * Read all of a File into a String.
   *
   * @param file the File to read.
   * @return a String containing the contents of the File
   * @exception IOException If there is a problem during reading.
   */
  public static String readAll(final File file) throws IOException {
    try (InputStream input = new FileInputStream(file)) {
      return readAll(input);
    }
  }


  /**
   * Read all of an input stream into a String.
   *
   * @param input input stream being read.
   * @return a String containing the contents of the input stream.
   * @exception IOException If there is a problem during reading.
   */
  public static String readAll(final InputStream input) throws IOException {
    return readAll(new InputStreamReader(input));
  }


  /**
   * Read all of an input stream into a String with the specified character encoding.
   *
   * @param input input stream being read.
   * @param encoding the character encoding string.
   * @return a String containing the contents of the input stream.
   * @exception IOException If there is a problem during reading.
   */
  public static String readAll(final InputStream input, final String encoding) throws IOException {
    return readAll(new InputStreamReader(input, encoding));
  }

  /**
   * Read all of a Reader into a String.
   *
   * @param input Reader being read.
   * @return a String containing the contents of the input stream.
   * @exception IOException If there is a problem during reading.
   */
  public static String readAll(final Reader input) throws IOException {
    final char[] b = new char[BUFFER_LENGTH];
    try (StringWriter str = new StringWriter(BUFFER_LENGTH)) {
      while (true) {
        final int length = input.read(b);
        if (length == EOF) {
          break;
        } else if (length == 0) {
          throw new RuntimeException();
        } else {
          str.write(b, 0, length);
        }
      }
      return str.toString();
    }
  }
}
