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
package com.rtg.ngs;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

import com.rtg.index.hash.ngs.OutputProcessor;
import com.rtg.launcher.HashingRegion;
import com.rtg.reader.NamesInterface;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.io.FileUtils;



/**
 * Creates a Top equals output processor, this will produce 'N' numbers of top results.
 *
 */
public class DefaultOutputProcessorSynch implements OutputProcessor {

  /** Name for temporary files */
  public static final String THREAD_FILE_NAME = "mapSubOutput_";

  private final NgsParams mParams;
  private final NamesInterface mNames;
  private final List<DefaultOutputProcessor> mChildren;

  /**
   * Create a new Top equal output processor
   * @param params ngs parameters
   * @throws IOException if an I/O Error occurs
   */
  public DefaultOutputProcessorSynch(final NgsParams params) throws IOException {
    mParams = params;
    if (mParams.searchParams() == null) {
      mNames = null;
    } else {
      mNames = mParams.searchParams().reader().names();
    }
    mChildren = new ArrayList<>();
  }

  @Override
  public void process(final long templateId, final String frame, final int readId, final int start, final int score, final int scoreIndel) {
    throw new UnsupportedOperationException("Should be called on thread local clones");
  }

  @Override
  public void close() throws IOException {
    for (final DefaultOutputProcessor child : mChildren) {
      child.close();
    }
  }

  private static class DefaultOutputProcessorComparator implements Comparator<DefaultOutputProcessor>, Serializable {
    @Override
      public int compare(final DefaultOutputProcessor first, final DefaultOutputProcessor other) {
      return first.getRegion().compareTo(other.getRegion());
    }
  }

  @Override
  public void finish() throws IOException {

    Collections.sort(mChildren, new DefaultOutputProcessorComparator());
    final File[] outputFiles = new File[mChildren.size()];
    for (int i = 0; i < outputFiles.length; ++i) {
      outputFiles[i] = mChildren.get(i).getFile();
    }
    if (mParams.outputParams().mergeMatchResults()) {
      try (OutputStream stream = mParams.outputParams().outStream()) {
        appendAllFiles(stream, outputFiles, true);
      }
      deleteTemp(mParams.outputParams().directory(), mChildren.size());
    }
  }

  static void appendAllFiles(final OutputStream dest, final File[] files, final boolean tempZip) throws IOException {
    boolean writeHeader = true;
    for (final File file : files) {
      writeHeader = !appendFile(dest, file, tempZip, writeHeader) && writeHeader;
    }
  }

  static void deleteTemp(final File inputDir, final int size) {
    for (int i = 0; i < size; ++i) {
      final File file = new File(inputDir, THREAD_FILE_NAME + i);
      if (!file.delete()) {
        Diagnostic.warning("Failed to delete temporary file: " + file.getPath());
      }
    }
  }

  private static boolean appendFile(final OutputStream out, final File in, final boolean gzipIn, final boolean copyHeader) throws IOException {
    //System.err.println("in=" + in.getAbsolutePath());
    boolean ret = false; //wrote out data
    if (in.length() > 0) {
      boolean skipLine = !copyHeader;
      try (InputStream inStream = gzipIn ? FileUtils.createGzipInputStream(in, true) : FileUtils.createFileInputStream(in, true)) {
        final byte[] buf = new byte[4096];
        int len;
        while ((len = inStream.read(buf, 0, buf.length)) > 0) {
          if (skipLine) {
            for (int i = 0; i < len; ++i) {
              if (buf[i] == (byte) '\n') {
                out.write(buf, i + 1, len - i - 1);
                skipLine = false;
                ret = true;
                break;
              }
            }
          } else {
            out.write(buf, 0, len);
            ret = true;
          }
        }

      }
    }
    return ret;
  }

  private static synchronized void createDir(final File dir) throws IOException {
    if (!dir.exists() && !dir.mkdirs()) {
      throw new IOException("unable to create directory: " + dir);
    }
  }

  @Override
  public synchronized OutputProcessor threadClone(final HashingRegion region) throws IOException {
    final int current = mChildren.size();
    final File dir = mParams.outputParams().directory();
    createDir(dir);
    final File outputFile = new File(dir, THREAD_FILE_NAME + current);
    //System.err.println("out=" + outputFile.getAbsolutePath());
    final OutputStream stream = FileUtils.createOutputStream(outputFile, true);
    final DefaultOutputProcessor child = new DefaultOutputProcessor(stream, mNames, region, outputFile);
    mChildren.add(child);
    if (region != HashingRegion.NONE) {
      return new ClippedOutputProcessor(child, region);
    } else {
      return child;
    }
  }

  @Override
  public void threadFinish() throws IOException {
    try {
      finish();
    } finally {
      close();
    }
  }
}
