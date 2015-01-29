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


import java.io.BufferedInputStream;
import java.io.BufferedReader;
import java.io.EOFException;
import java.io.File;
import java.io.FileDescriptor;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.io.Reader;
import java.util.Locale;
import java.util.zip.GZIPInputStream;

import com.rtg.util.Resources;
import com.rtg.util.Utils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.ErrorType;
import com.rtg.util.diagnostic.NoTalkbackSlimException;
import com.rtg.util.diagnostic.OneShotTimer;
import com.rtg.util.gzip.GzipUtils;
import com.rtg.util.io.bzip2.CBZip2InputStream;

/**
 * Utility functions for manipulating files that are not provided in the File
 * class.
 *
 */
public final class FileUtils {

  /** The suitable size for a buffer for file i/o */
  public static final int BUFFER_SIZE = 8192;

  // Buffer size, close to the disk
  /**
   * Buffer size to be used for buffered streams.
   */
  public static final int BUFFERED_STREAM_SIZE = 1024 * 1024; // 1MB upsets disks

  /**
   * Name of the log file.
   */
  public static final String LOG_SUFFIX = "log";

  /**
   * Name of the progress file.
   */
  public static final String PROGRESS_SUFFIX = "progress";

  /**
   * Name of the suffix used on compressed output files
   */
  public static final String GZ_SUFFIX = ".gz";

  /**
   * Name of the bzip2 suffix
   */
  public static final String BZ2_SUFFIX = ".bz2";

  private FileUtils() {
  }

  /**
   * Patches a string URL so that it works for UNC files. This is a work-around
   * for a bug in Java's handling of UNC files. URL's created from a string such
   * as:
   *
   * <pre><code>
   * new URL(&quot;file://10.65.1.7/path/to/file&quot;);
   * </code></pre>
   *
   * cannot be read, even though that string representation is what is returned
   * by the functioning URL returned by:
   *
   * <pre><code>
   *    new File(&quot;\\10.65.1.7\path\to\file&quot;).toURL();
   * </code></pre>
   *
   * It turns out that the following URL would work:
   *
   * <pre><code>
   * new URL(&quot;file://\\10.65.1.7/path/to/file&quot;);
   * </code></pre>
   *
   * so that is what this method returns, if the URL starts with "file://".
   *
   * @param url a <code>String</code> representation of a URL
   * @return a <code>String</code> representation that will function for UNC
   * URLs
   */
  public static String patchURLForUNC(final String url) {
    return ((url.length() > 8) && (url.startsWith("file://")) && (url.charAt(7) != '/') && (url
        .charAt(7) != '\\')) ? "file://\\\\" + url.substring(7) : url;
  }

  /**
   * Checks to see if the given file name corresponds to a file that exists and
   * can be read.
   *
   * @param fileName a file name to test.
   * @return true iff fileName exists and can be read.
   */
  public static boolean checkFile(final String fileName) {
    return checkFile(new File(fileName));
  }

  /**
   * Checks to see if the given file exists, is a file, and can be read.
   *
   * @param f a File to test.
   * @return true iff f exists, is a file (not a directory) and can be read.
   */
  public static boolean checkFile(final File f) {
    return f.exists() && f.isFile() && f.canRead();
  }

  /**
   * Checks if the given file exists, can be read and is not empty.
   *
   * @param fileName the filename of the File to test
   * @return true iff f is nonempty, exists and can be read.
   */
  public static boolean isNonEmpty(final String fileName) {
    return isNonEmpty(new File(fileName));
  }

  /**
   * Checks if the given file exists, can be read and is not empty.
   *
   * @param f the File to test
   * @return true iff f is nonempty, exists and can be read.
   */
  public static boolean isNonEmpty(final File f) {
    return checkFile(f) && f.length() > 0;
  }

  /**
   * Checks to see if the given directory name corresponds to a directory that
   * exists and can be read.
   *
   * @param dirName a file name to test.
   * @return true iff dirName exists, is a directory and can be read.
   */
  public static boolean checkDir(final String dirName) {
    return checkDir(new File(dirName));
  }

  /**
   * Checks to see if the given file exists and is a directory.
   *
   * @param f a File to test.
   * @return true iff f exists and is a directory.
   */
  public static boolean checkDir(final File f) {
    return f.exists() && f.isDirectory();
  }

  /**
   * Concatenates multiple files together, no frills.
   *
   * @param filename destination for cat
   * @param deleteIntermediate delete intermediate files as you are done with them
   * @param inputFiles files to merge, in order
   * @throws IOException if an I/O error occurs
   */
  public static void catInSync(final File filename, boolean deleteIntermediate, final File... inputFiles) throws IOException {
    final byte[] buff = new byte[100 * 1024]; //100k for our buffer, as we sync on the FileOutputStream.
    final OneShotTimer timer = new OneShotTimer("catInSync");
    try (FileOutputStream destination = new FileOutputStream(filename)) {
      for (File inputFile : inputFiles) {
        long t0 = System.nanoTime();
        final long start = System.nanoTime();
        final long length = inputFile.length();
        Diagnostic.developerLog("start catInSync file=" + inputFile.getAbsolutePath() + " bytes=" + length);
        if (length > 0) {
          try (InputStream input = new AsynchInputStream(FileUtils.createFileInputStream(inputFile, true))) {
            int len;
            //          double last = System.nanoTime() / 1000000000.0; // current time since start
            while ((len = input.read(buff)) > 0) {
              destination.write(buff, 0, len);
              if (System.nanoTime() - t0 > 5000000000L) {
                // sync every 5 seconds
                destination.getFD().sync(); // we are doing this at the moment so that the disk is completely flushed to disk between buffer
                // writes. We are really guessing, but it appears to have made a differences in manky drives.
                t0 = System.nanoTime();
              }
              //            final double cur = System.nanoTime() / 1000000000.0; // current time since start
              // log if less than 1MB / second
              //            if (len / (cur - last) < 1000000) {
              //              Diagnostic.developerLog("[ " + last + " ] last time " + len);
              //              Diagnostic.developerLog("[ " + cur + " ] finished read/write " + len + " speed= " + (len / (cur - last)));
              //            }
              //           last = cur;
            }
          }
        }
        final long diff = System.nanoTime() - start;
        Diagnostic.developerLog("end catInSync file=" + inputFile.getAbsolutePath() + " bytes=" + length
          + " time=" + (diff / 1000000) + "ms"
          + " bytes/sec=" + Utils.realFormat(length * 1.0e9 / diff, 2));
        if (deleteIntermediate) {
          if (!inputFile.delete()) {
            Diagnostic.userLog("Failed to delete intermediate file: " + inputFile.getPath());
          }
        }
      }
    }
    timer.stopLog();
  }



  /**
   * Creates an empty directory in the default temporary-file directory, using
   * the given prefix and suffix to generate its name. Invoking this method is
   * equivalent to invoking <code>createTempDir(prefix, suffix, null)</code>.
   *
   * @param prefix The prefix string to be used in generating the file's name;
   * must be at least three characters long
   * @param suffix The suffix string to be used in generating the file's name;
   * may be null, in which case the suffix "<code>.tmp</code>" will be used
   * @return An abstract pathname denoting a newly-created empty directory
   * @exception IOException If a directory could not be created
   * @exception IllegalArgumentException If the prefix argument contains fewer
   * than three characters
   * @exception SecurityException If a security manager exists and its
   * <code>SecurityManager.checkWrite(java.lang.String)</code> method does not allow a file
   * to be created
   */
  public static File createTempDir(final String prefix, final String suffix) throws IOException {
    return createTempDir(prefix, suffix, null);
  }

  /**
   * Creates a new empty directory in the specified directory, using the given
   * prefix and suffix strings to generate its name.
   *
   * @param prefix The prefix string to be used in generating the file's name;
   * must be at least three characters long
   * @param suffix The suffix string to be used in generating the file's name;
   * may be null, in which case the suffix "<code>.tmp</code>" will be used
   * @param directory The directory in which the file is to be created, or null
   * if the default temporary-file directory is to be used
   * @return An abstract pathname denoting a newly-created empty directory
   * @exception IOException If a directory could not be created
   * @exception IllegalArgumentException If the prefix argument contains fewer
   * than three characters
   * @exception SecurityException If a security manager exists and its
   * <code>SecurityManager.checkWrite(java.lang.String</code>) method does not allow a file
   * to be created
   */
  public static File createTempDir(final String prefix, final String suffix, final File directory)
      throws IOException {
    final File dir = File.createTempFile(prefix, suffix, directory);
    if (!dir.delete()) {
      throw new IOException("Could not create directory: " + dir.toString());
    }
    if (!dir.mkdir()) {
      throw new IOException("Could not create directory: " + dir.toString());
    }
    return dir;
  }

  /**
   * Creates a file copying the string <code>content</code> to it.
   *
   * @param content a <code>String</code>
   * @param file a <code>File</code> to write to
   * @return a <code>File</code> containing the string content
   * @exception IOException if an error occurs.
   * @exception NullPointerException if the content is null
   */
  public static File stringToFile(final String content, final File file) throws IOException {
    if (content == null) {
      throw new NullPointerException("null string given");
    }
    return byteArrayToFile(content.getBytes(), file);
  }

  /**
   * Appends the contents String to a file.
   *
   * @param file the file to save to.
   * @param contents String to append to file.
   * @exception IOException if there is an error opening or writing to file.
   */
  public static void appendToFile(final File file, final String contents) throws IOException {
    try (FileOutputStream fos = new FileOutputStream(file, true)) {
      fos.write(contents.getBytes());
    }
  }

  /**
   * Creates a file copying the byte array <code>content</code> to it.
   *
   * @param content a array of <code>byte</code>s
   * @param file a <code>File</code> to write to
   * @return a <code>File</code> containing the byte array content
   * @exception IOException if an error occurs.
   * @exception NullPointerException if the content is null
   */
  public static File byteArrayToFile(final byte[] content, final File file) throws IOException {
    if (content == null) {
      throw new NullPointerException("null array of bytes given");
    }
    if (file == null) {
      throw new NullPointerException("null file given");
    }
    try (FileOutputStream out = new FileOutputStream(file)) {
      out.write(content);
    }
    return file;
  }

  /**
   * Creates an output stream for a file, determining from the filename whether to compress.
   *
   * @param file the output <code>File</code>
   * @return an <code>OutputStream</code> value
   * @exception IOException if an error occurs.
   */
  public static OutputStream createOutputStream(File file) throws IOException {
    return createOutputStream(file, FileUtils.isGzipFilename(file), false);
  }

  /**
   * Creates an output stream for a file, wrapping in a
   * <code>BufferedOutputStream</code> stream and optionally a
   * <code>GZIPOutputStream</code>.
   *
   * @param file the output <code>File</code>
   * @param zip if true, the output will be gzip compressed.
   * @return an <code>OutputStream</code> value
   * @exception IOException if an error occurs.
   */
  public static OutputStream createOutputStream(File file, boolean zip) throws IOException {
    return createOutputStream(file, zip, false);
  }

  /**
   * Creates an output stream for a file, wrapping in a
   * <code>BufferedOutputStream</code> OutputStream and optionally a
   * <code>GZIPOutputStream</code>.
   *
   * @param file the output <code>File</code>
   * @param zip if true, the output will be gzip compressed.
   * @param append if true, the output will be appended to the file.
   * @return an <code>OutputStream</code> value
   * @exception IOException if an error occurs.
   */
  public static OutputStream createOutputStream(File file, boolean zip, boolean append) throws IOException {
    return createOutputStream(file, zip, append, false);
  }

  /**
   * Creates an output stream for a file, wrapping in a
   * <code>BufferedOutputStream</code> OutputStream and optionally a
   * <code>GZIPOutputStream</code>.
   *
   * @param file the output <code>File</code>
   * @param zip if true, the output will be gzip compressed.
   * @param append if true, the output will be appended to the file.
   * @param terminate if true (and block compressed), terminate the file
   * @return an <code>OutputStream</code> value
   * @exception IOException if an error occurs.
   */
  public static OutputStream createOutputStream(File file, boolean zip, boolean append, boolean terminate) throws IOException {
    Diagnostic.developerLog("FileUtils.outputStream " + file.getAbsolutePath() + " " + zip + " " + append);
    OutputStream outStream = new FileOutputStream(file, append);
    outStream = new BufferedOutputStreamFix(outStream, BUFFERED_STREAM_SIZE);
    if (zip) {
      outStream = new GzipAsynchOutputStream(outStream, terminate);
    }
    return outStream;
  }


  /**
   * Creates a <code>Reader</code> for a file automatically
   * choosing between basic file input, gzip input and bzip2 input based on
   * file extension.
   *
   * @param file the input <code>File</code>
   * @param async if true, the input will be asynchronous (not to be used with picard)
   * @return a <code>Reader</code> value
   * @exception IOException if an error occurs.
   */
  public static Reader createReader(File file, boolean async) throws IOException {
    final InputStream is = createInputStream(file, async);
    return new InputStreamReader(is);
  }

  /**
   * Creates a <code>BufferedInputStream</code> for a file automatically
   * choosing between basic file input, gzip input and bzip2 input based on
   * file extension.
   *
   * @param file the input <code>File</code>
   * @param async if true, the input will be asynchronous (not to be used with picard)
   * @return an <code>BufferedInputStream</code> value
   * @exception IOException if an error occurs.
   */
  public static BufferedInputStream createInputStream(File file, boolean async) throws IOException {
    if (FileUtils.isGzipFilename(file)) {
      return createGzipInputStream(file, async);
    } else if (FileUtils.isBzip2Filename(file)) {
      return createBzip2InputStream(file, async);
    } else {
      return createFileInputStream(file, async);
    }
  }

  /**
   * Creates a <code>BufferedInputStream</code> for a file.
   *
   * @param file the input <code>File</code>
   * @param async if true, the input will be asynchronous (not to be used with picard)
   * @return an <code>BufferedInputStream</code> value
   * @exception IOException if an error occurs.
   */
  public static BufferedInputStream createFileInputStream(File file, boolean async) throws IOException {
    final InputStream inStream;
    if (async) {
      inStream = new AsynchInputStream(new FileInputStream(file));
    } else {
      inStream = new FileInputStream(file);
    }
    return new BufferedInputStream(inStream, BUFFERED_STREAM_SIZE);
  }

  /**
   * Creates a <code>BufferedInputStream</code> for a file,
   * utilizing a <code>GZIPInputStream</code>.
   *
   * @param file the input <code>File</code>
   * @param async if true, the input will be asynchronous (not to be used with picard)
   * @return an <code>BufferedInputStream</code> value
   * @throws IOException if an error occurs.
   */
  public static BufferedInputStream createGzipInputStream(File file, boolean async) throws IOException {
    final InputStream inStream;
    if (async) {
      inStream = new GzipAsynchInputStream(file);
    } else {
      inStream = GzipUtils.createGzipInputStream(new BufferedInputStream(new FileInputStream(file)));
    }
    return new BufferedInputStream(inStream);
  }

  /**
   * Creates a <code>BufferedInputStream</code> for a file,
   * utilizing a <code>CBZip2InputStream</code>.
   *
   * @param file the input <code>File</code>
   * @param async if true, the input will be asynchronous (not to be used with picard)
   * @return an <code>BufferedInputStream</code> value
   * @throws IOException if an error occurs.
   */
  public static BufferedInputStream createBzip2InputStream(File file, boolean async) throws IOException {
    final InputStream inStream;
    if (async) {
      inStream = new AsynchInputStream(new CBZip2InputStream(new BufferedInputStream(new FileInputStream(file), FileUtils.BUFFERED_STREAM_SIZE)));
    } else {
      inStream = new CBZip2InputStream(new BufferedInputStream(new FileInputStream(file), FileUtils.BUFFERED_STREAM_SIZE));
    }
    return new BufferedInputStream(inStream);
  }


  /**
   * Test if the supplied directory is valid for the purposes of writing
   * a SLIM result.  This means the directory either exists or
   * does not exist and can be successfully created
   * (which this method will do).
   *
   * @param directory directory to test
   */
  public static void ensureOutputDirectory(final File directory) {
    if (!directory.exists()) {
      if (!directory.mkdirs()) {
        throw new NoTalkbackSlimException(ErrorType.DIRECTORY_NOT_CREATED, directory.getPath());
      }
    } else {
      if (!directory.isDirectory()) {
        throw new NoTalkbackSlimException(ErrorType.NOT_A_DIRECTORY, directory.getPath());
      }
    }
  }

  /**
   * Read the contents of a file and turn it into a string.
   * @param file to be read.
   * @return the contents of the file.
   * @throws IOException if an I/O error occurs.
   */
  public static String fileToString(final String file) throws IOException {
    final String ret;
    try (Reader reader = new FileReader(file)) {
      ret = readerToString(reader);
    }
    return ret;
  }

  /**
   * Read the contents of a file and turn it into a string.
   * @param file to be read.
   * @return the contents of the file.
   * @throws IOException if an I/O error occurs.
   */
  public static String fileToString(final File file) throws IOException {
    final String ret;
    try (Reader reader = new FileReader(file)) {
      ret = readerToString(reader);
    }
    return ret;
  }

  /**
   * Read the contents of a reader and turn it into a string.
   * @param fileReader reader.
   * @return the contents of the reader.
   * @throws IOException if an I/O error occurs.
   */
  public static String readerToString(final Reader fileReader) throws IOException {
    final StringBuilder sb = new StringBuilder();
    try (BufferedReader br = new BufferedReader(fileReader)) {
      final char[] buffer = FileUtils.makeBuffer();
      final int eof = -1;
      for (int len = br.read(buffer); len > eof; len = br.read(buffer)) {
        for (int i = 0; i < len; i++) {
          sb.append(buffer[i]);
        }
      }
    }
    return sb.toString();
  }

  /**
   * @return a buffer of suitable size for file i/o.
   */
  public static char[] makeBuffer() {
    return new char[FileUtils.BUFFER_SIZE];
  }

  /**
   * Gets the contents of the stream as a string.
   *
   * @param stream an <code>InputStream</code>
   * @return a String containing the contents of the stream.
   * @exception IOException if an error occurs.
   * @exception NullPointerException if the stream is null
   */
  public static String streamToString(final InputStream stream) throws IOException {
    if (stream == null) {
      throw new NullPointerException("null stream given");
    }
    final StringBuilder out = new StringBuilder();
    final byte[] b = FileUtils.makeByteBuffer();
    int len = stream.read(b);
    final int eof = -1;
    while (len != eof) {
      for (int i = 0; i < len; i++) {
        out.append((char) b[i]); //this does not do careful conversion but good enough for testing use.
      }
      len = stream.read(b);
    }
    return out.toString();
  }

  /**
   * @return a buffer of suitable size for file i/o.
   */
  static byte[] makeByteBuffer() {
    return new byte[FileUtils.BUFFER_SIZE];
  }

  /**
   * Get current process id.
   *
   * @return process id.
   */
  public static int getProcessId() {
    return Integer.parseInt(java.lang.management.ManagementFactory.getRuntimeMXBean().getName().split("@")[0]);
  }

  private static OutputStream javaGetStdoutAsOutputStream() {
    return new FileOutputStream(FileDescriptor.out) {
      @Override
      public void close() throws IOException {
        // Bad stuff happens if you really close stdout
        flush();
      }
    };
  }

  /**
   * Return the standard output stream as an ordinary output stream rather than
   * a <code>PrintStream</code>.
   *
   * @return standard output stream
   */
  public static OutputStream getStdoutAsOutputStream() {
    return javaGetStdoutAsOutputStream();
  }

  private static OutputStream javaGetStderrAsOutputStream() {
    return new FileOutputStream(FileDescriptor.err) {
      @Override
      public void close() throws IOException {
        // Bad stuff happens if you really close stderr
        flush();
      }
    };
  }

  /**
   * Return the standard error stream as an ordinary output stream rather than
   * a <code>PrintStream</code>.
   *
   * @return standard error stream
   */
  public static OutputStream getStderrAsOutputStream() {
    return javaGetStderrAsOutputStream();
  }

  /**
   * Utility function for deleting the files in a given directory and the directory itself.
   * (Files only, non-recursive).
   * @param file the directory to delete the files from.
   * @return true if all the files were successfully deleted.
   */
  public static boolean deleteFiles(final File file) {
    boolean ok = true;
    if (file != null) {
      if (file.isDirectory()) {
        final File[] files = file.listFiles();
        if (files != null) {
          for (final File f : files) {
            if (f.isFile()) {
              ok &= f.delete();
            }
          }
        } else {
          ok = false;
        }
        ok &= file.delete();
      }
    }
    return ok;
  }

  /**
   * A wrapper for stream.skip that is guaranteed to skip at least 1 byte when called
   * @param stream the stream to skip along
   * @param count how far to skip
   * @return distance skipped
   * @throws IOException if an IOException occurs
   */
  public static long streamSkip(InputStream stream, long count) throws IOException {
    long skipped = stream.skip(count);
    if (skipped < 1) {
      // Skip has failed lets attempt to use read to fix things up.
      final byte[] temp = new byte[1024];
      skipped = stream.read(temp, 0, Math.min((int) count, 1024));
    }
    return skipped;
  }

  /**
   * A wrapper for stream.skip that is guaranteed to skip the specified amount
   * or throw an IOException if that is not possible
   * @param is the stream to skip along
   * @param amount how far to skip
   * @throws IOException if an IOException occurs
   */
  public static void skip(InputStream is, long amount) throws IOException {
    long remaining = amount;
    long read;
    while (remaining > 0 && (read = FileUtils.streamSkip(is, remaining)) > 0) {
      remaining -= read;
    }
    if (remaining > 0) {
      throw new EOFException();
    }
  }

  /**
   * Loads an integer from a file
   * @param file the file containing the integers
   * @param index the index of the integer to get
   * @return the integer at <code>index</code> in <code>file</code>
   * @throws IOException if an IOException occurs
   */
  public static int getIntFromFile(File file, int index) throws IOException {
    final byte[] bytes = new byte[4];
    try (BufferedInputStream bis = new BufferedInputStream(new FileInputStream(file))) {
      FileUtils.streamSkip(bis, (long) index * 4);
      int remaining = 4;
      int read;
      while (remaining > 0 && (read = bis.read(bytes, bytes.length - remaining, remaining)) > 0) {
        remaining -= read;
      }
      if (remaining > 0) {
        throw new EOFException();
      }
    }
    return ByteArrayIOUtils.bytesToIntBigEndian(bytes, 0);
  }

  /**
   * @param file file to check
   * @return true if has <code>BZIP2</code> file extension
   */
  public static boolean isBzip2Filename(File file) {
    return isBzip2Filename(file.getName());
  }

  /**
   * @param file file to check
   * @return true if has <code>BZIP2</code> file extension
   */
  public static boolean isBzip2Filename(String file) {
    return file.toLowerCase(Locale.getDefault()).endsWith(BZ2_SUFFIX);
  }

  /**
   * @param file file to check
   * @return true if has GZIP file extension
   */
  public static boolean isGzipFilename(File file) {
    return isGzipFilename(file.getName());
  }

  /**
   * @param file file to check
   * @return true if has GZIP file extension
   */
  public static boolean isGzipFilename(String file) {
    return file.toLowerCase(Locale.getDefault()).endsWith(GZ_SUFFIX);
  }

  /**
   * Checks first 2 bytes for presence of GZIP identifier.
   * Stream is then reset to position it was in before being passed to this method.
   * @param in stream to check
   * @return true if contains GZIP identifier
   * @throws IOException if an IO error occurs
   */
  public static boolean isGzipFile(InputStream in) throws IOException {
    if (!in.markSupported()) {
      throw new IllegalArgumentException("Require a resetable stream");
    }
    in.mark(2);
    try {
      final byte[] b = new byte[2];
      IOUtils.readFully(in, b, 0, 2);
      final int magic = (((int) b[1] & 0xff) << 8) | ((int) b[0] & 0xff);
      return magic == GZIPInputStream.GZIP_MAGIC;
    } finally {
      in.reset();
    }
  }

  /**
   * Write the contents of <code>in</code> into <code>out</code>
   * @param in input stream
   * @param out output stream
   * @param bufSize size of the buffer to use
   * @throws IOException if an IO error occurs
   */
  public static void streamToStream(InputStream in, OutputStream out, int bufSize) throws IOException {
    final byte[] buf = new byte[bufSize];
    int len;
    while ((len = in.read(buf)) != -1) {
      out.write(buf, 0, len);
    }
  }

  /**
   * Copy a resource to a File
   * @param resource source resource
   * @param outFile destination file
   * @throws IOException if an IO error occurs
   */
  public static void copyResource(String resource, File outFile) throws IOException {
    try (FileOutputStream fos = new FileOutputStream(outFile)) {
      try (InputStream is = Resources.getResourceAsStream(FileUtils.class, resource)) {
        streamToStream(is, fos, FileUtils.BUFFER_SIZE);
      }
    }
  }

  /**
   * Converts a file into one with an appropriate extension, depending on whether the file should be gzipped.
   * @param gzip true if the output file is destined to be GZIP compressed.
   * @param file the input file
   * @return the appropriately adjusted file
   */
  public static File getZippedFileName(final boolean gzip, final File file) {
    if (gzip == isGzipFilename(file)) {    // All is good
      return file;
    } else if (gzip) {    // Need to add suffix
      return new File(file.getAbsolutePath() + GZ_SUFFIX);
    } else {    // Need to remove suffix
      final String path = file.getAbsolutePath();
      return new File(path.substring(0, path.length() - GZ_SUFFIX.length()));
    }
  }

  /**
   * Creates an output stream for a file, wrapping in a
   * <code>BufferedOutputStream</code> OutputStream and optionally a
   * <code>GZIPOutputStream</code>. Additionally all data being written to the
   * file will also be written to the other stream specified.
   *
   * @param file the output <code>File</code>
   * @param otherStream stream to copy output to
   * @param zip if true, the output will be gzip compressed.
   * @param append if true, the output will be appended to the file.
   * @param terminate if true (and block compressed), terminate the file
   * @return an <code>OutputStream</code> value
   * @exception IOException if an error occurs.
   */
  public static OutputStream createTeedOutputStream(File file, OutputStream otherStream, boolean zip, boolean append, boolean terminate) throws IOException {
    Diagnostic.developerLog("FileUtils.outputStream " + file.getAbsolutePath() + " " + zip + " " + append);
    OutputStream outStream = new FileOutputStream(file, append);
    outStream = new TeeOutputStream(new BufferedOutputStreamFix(outStream, BUFFERED_STREAM_SIZE), otherStream);
    if (zip) {
      outStream = new GzipAsynchOutputStream(outStream, terminate);
    }
    return outStream;
  }

  private static class TeeOutputStream extends OutputStream {
    private final OutputStream mOut1;
    private final OutputStream mOut2;

    public TeeOutputStream(OutputStream out1, OutputStream out2) {
      mOut1 = out1;
      mOut2 = out2;
    }

    @Override
    @SuppressWarnings("try")
    public void close() throws IOException {
      try (OutputStream ignored = mOut2) {
        mOut1.close();
      }
    }

    @Override
    public void flush() throws IOException {
      mOut1.flush();
      mOut2.flush();
    }

    @Override
    public void write(int b) throws IOException {
      mOut1.write(b);
      mOut2.write(b);
    }

    @Override
    public void write(byte[] b) throws IOException {
      mOut1.write(b);
      mOut2.write(b);
    }

    @Override
    public void write(byte[] b, int off, int len) throws IOException {
      mOut1.write(b, off, len);
      mOut2.write(b, off, len);
    }

  }

}
