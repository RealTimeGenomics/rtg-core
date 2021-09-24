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
