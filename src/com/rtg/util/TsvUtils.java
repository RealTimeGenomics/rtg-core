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

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;

import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.OneShotTimer;
import com.rtg.util.io.FileUtils;

/**
 * Functions to manipulate tab delimited files.
 *
 */
public final class TsvUtils {

  private TsvUtils() {

  }

  /**
   * Cats multiple tab separated value files together
   * @param destination destination for cat
   * @param inputFiles files to merge, in order
   * @throws IOException if an I/O error occurs
   */
  public static void tsvCat(OutputStream destination, File... inputFiles) throws IOException {
    tsvCat(inputFiles.length > 0 && FileUtils.isGzipFilename(inputFiles[0]), destination, inputFiles);
  }

  /**
   * Cats multiple tab separated value files together
   * @param gzipped true if the input files are gzipped
   * @param destination destination for cat
   * @param inputFiles files to merge, in order
   * @throws IOException if an I/O error occurs
   */
  public static void tsvCat(boolean gzipped, OutputStream destination, File... inputFiles) throws IOException {
    final byte[] buff = new byte[1024 * 1024];
    final OneShotTimer timer = new OneShotTimer("tsvCat");
    for (int i = 0; i < inputFiles.length; i++) {
      final long t0 = System.nanoTime();
      boolean dropHeader = i > 0;
      boolean scanNewLine = false;
      final File current = inputFiles[i];
      final long length = current.length();
      if (length > 0) {
        try (InputStream input = gzipped ? FileUtils.createGzipInputStream(current, true) : FileUtils.createFileInputStream(current, true)) {
          int len;
          while ((len = input.read(buff)) > 0) {
            int currentPos = 0;
            if (dropHeader) {
              for (; currentPos < len; currentPos++) {
                final char c = (char) buff[currentPos];
                if (scanNewLine && c == '\n') {
                  scanNewLine = false;
                } else if (!scanNewLine) {
                  if (c == '#') {
                    scanNewLine = true;
                  } else {
                    dropHeader = false;
                    break;
                  }
                }
              }
            }
            //if dropHeader == true here, then len should == currentPos
            assert dropHeader == (len == currentPos);
            destination.write(buff, currentPos, len - currentPos);
          }
        }
      }
      final long diff = System.nanoTime() - t0;
      Diagnostic.developerLog("tsv concat file=" + current.getAbsolutePath() + " bytes=" + length
                              + " time=" + (diff / 1000000) + "ms"
                              + " bytes/sec=" + Utils.realFormat(length * 1.0e9 / diff, 2));
    }
    timer.stopLog();
  }

}
