/*
 * Copyright (c) 2018. Real Time Genomics Limited.
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the
 *    distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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
   * Cats multiple tab separated value files together, keeping the header from only the first file.
   * @param gzipped true if the input files are gzipped
   * @param destination destination for cat
   * @param inputFiles files to merge, in order
   * @throws IOException if an I/O error occurs
   */
  public static void tsvCat(boolean gzipped, OutputStream destination, File... inputFiles) throws IOException {
    final byte[] buff = new byte[1024 * 1024];
    final OneShotTimer timer = new OneShotTimer("tsvCat");
    for (int i = 0; i < inputFiles.length; ++i) {
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
              for (; currentPos < len; ++currentPos) {
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
