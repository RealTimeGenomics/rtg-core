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

package com.rtg.variant.sv;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.Arrays;

import com.rtg.util.integrity.Exam;
import com.rtg.util.integrity.IntegralAbstract;
import com.rtg.util.io.FileUtils;

/**
 * Maintain information about corrections to be applied to coverage counts.
 */
//TODO allow multiple sequences and use tabix indexes.
public class Corrections extends IntegralAbstract {

  private double[] mCorrections;

  /**
   * @param file with corrections in tab separated format (sequence name, position - 0 based, correction - float &ge; 0.0). Likely to be
   * <code>bgzipped</code>.
   * @throws IOException whenever.
   */
  Corrections(final File file) throws IOException {
    this(FileUtils.createGzipInputStream(file, true), 1000);
  }

  /**
   * For use only in testing.
   * @param inStr input stream with corrections.
   * @param initialLength initial length of array.
   * @throws IOException whenever.
   */
  Corrections(final InputStream inStr, int initialLength) throws IOException {
    mCorrections = new double[initialLength];
    final BufferedReader in = new BufferedReader(new InputStreamReader(inStr));
    String seq = null;
    while (true) {
      final String line = in.readLine();
      if (line == null) {
        break;
      }
      if (line.length() == 0) {
        continue;
      }
      if (line.charAt(0) == '#') {
        continue;
      }
      final String[] split = line.split("\\s+");
      if (split.length != 3) {
        throw new RuntimeException(line);
      }
      if (seq == null) {
        seq = split[0];
      } else {
        if (!seq.equals(split[0])) {
          throw new RuntimeException(line);
        }
      }
      try {
        final int index = Integer.parseInt(split[1]);
        final double v = Double.parseDouble(split[2]);
        add(index, v);
      } catch (final Exception e) {
        throw new RuntimeException(line, e);
      }
    }
    in.close();
  }

  // for testing only.
  int length() {
    return mCorrections.length;
  }
  private void add(final int index, final double v) {
    //System.err.println("index=" + index + " length=" + mCorrections.length);
    if (index >= mCorrections.length) {
      final int nuLength = ((Math.max(index, mCorrections.length) + 1) * 3) / 2;
      mCorrections = Arrays.copyOf(mCorrections, nuLength);
      //System.err.println("resize length=" + mCorrections.length);
    }
    mCorrections[index] = v;
  }

  /**
   * Get the correction for current sequence.
   * @param index position for correction.
   * @return correction (&ge; 0.0).
   */
  double correction(final int index) {
    return mCorrections[index]; //TODO
  }

  @Override
  public boolean integrity() {
    Exam.assertNotNull(mCorrections);
    Exam.assertTrue(mCorrections.length >= 1000);
    return true;
  }
}
