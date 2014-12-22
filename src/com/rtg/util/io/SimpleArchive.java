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

import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;

/**
 * Simple archive format, doesn't support directories
 */
public final class SimpleArchive {

  private static final int VERSION = 1;
  private SimpleArchive() {
  }

  private static class FileHeader {
    final int mNameSize;
    final long mFileSize;

    public FileHeader(int nameSize, long fileSize) {
      this.mNameSize = nameSize;
      this.mFileSize = fileSize;
    }
  }

  /**
   * extracts archive
   * @param archive archive file
   * @param outputDir destination dir
   * @throws IOException if an IO error occurs
   */
  public static void unpackArchive(File archive, File outputDir) throws IOException {
    try (DataInputStream input = new DataInputStream(new FileInputStream(archive))) {
      unpackArchivePrivate(input, outputDir);
    }
  }

  /**
   * extracts archive
   * @param input stream from archive file
   * @param outputDir destination dir
   * @throws IOException if an IO error occurs
   */
  public static void unpackArchive(InputStream input, File outputDir) throws IOException {
    final DataInputStream dInput = new DataInputStream(input);
    unpackArchivePrivate(dInput, outputDir);
  }

  /**
   * extracts archive
   * @param input stream from archive file
   * @param outputDir destination dir
   * @throws IOException if an IO error occurs
   */
  private static void unpackArchivePrivate(DataInputStream input, File outputDir) throws IOException {
    if (!outputDir.exists()) {
      if (!outputDir.mkdirs()) {
        throw new IllegalArgumentException("Could not create directory: " + outputDir.getPath());
      }
    }
    input.readInt(); //version
    final byte[] buf = new byte[4096];
    final int numFiles = input.readInt();
    for (int i = 0; i < numFiles; i++) {
      final FileHeader header = new FileHeader(input.readInt(), input.readLong());
      final byte[] name = new byte[header.mNameSize];
      input.readFully(name);
      final File outFile = new File(outputDir, new String(name));
      try (FileOutputStream output = new FileOutputStream(outFile)) {
        long remaining = header.mFileSize;
        int toRead = buf.length < remaining ? buf.length : (int) remaining;
        int len;
        while (toRead > 0 && (len = input.read(buf, 0, toRead)) > 0) {
          output.write(buf, 0, len);
          remaining -= len;
          toRead = buf.length < remaining ? buf.length : (int) remaining;
        }
      }
    }
  }

  /**
   * Creates archive
   * @param dest file to contain other files
   * @param input input files
   * @throws IOException if an IO error occurs
   */
  public static void writeArchive(File dest, File... input) throws IOException {
    try (DataOutputStream output = new DataOutputStream(new FileOutputStream(dest))) {
      output.writeInt(VERSION);
      output.writeInt(input.length);
      final byte[] buf = new byte[4096];
      for (File f : input) {
        final byte[] nameBytes = f.getName().getBytes();
        final int nameSize = nameBytes.length;
        final FileHeader header = new FileHeader(nameSize, f.length());

        output.writeInt(header.mNameSize);
        output.writeLong(header.mFileSize);
        output.write(nameBytes);

        final FileInputStream fis = new FileInputStream(f);
        try {
          int len;
          while ((len = fis.read(buf)) > 0) {
            output.write(buf, 0, len);
          }
        } finally {
          fis.close();
        }
      }
    }
  }

  private static void printUsage() {
    System.out.println("Usage: ");
    System.out.println("SimpleArchive c archive files...");
    System.out.println("SimpleArchive d archive outputDir");
  }

  /**
   * Creates or unpacks archive
   * @param args arguments
   * @throws IOException if an io error occurs
   */
  public static void main(String[] args) throws IOException {
    if (args.length < 3) {
      printUsage();
      return;
    }
    if (args[0].equals("c")) {
      final File archive = new File(args[1]);
      final File[] input = new File[args.length - 2];
      for (int i = 0; i < input.length; i++) {
        input[i] = new File(args[i + 2]);
      }
      writeArchive(archive, input);
    } else if (args[0].equals("d")) {
      final File archive = new File(args[1]);
      final File output = new File(args[2]);
      unpackArchive(archive, output);
    } else {
     printUsage();
    }
  }
}
