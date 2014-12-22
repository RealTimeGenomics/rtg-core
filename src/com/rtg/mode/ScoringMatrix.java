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
package com.rtg.mode;

import static com.rtg.util.StringUtils.LS;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.Serializable;
import java.util.Arrays;

import com.rtg.util.integrity.Exam;
import com.rtg.util.integrity.Integrity;

/**
 * An abstract class to keep common components of  Blosum62.java and ScoringMatrix.java
 */
public abstract class ScoringMatrix implements Serializable, Integrity {

  protected int[][] mScores;

  /**
   * Get the integer value of the score.
   * @param i first ordinal value of an amino acid.
   * @param j second ordinal value of an amino acid.
   * @return score.
   */
  public int score(final int i, final int j) {
    return mScores[i][j];
  }

  @Override
  public String toString() {
    final StringBuilder sb = new StringBuilder();
    for (final int[] mScore : mScores) {
      sb.append(Arrays.toString(mScore)).append(LS);
    }
    return sb.toString();
  }

  @Override
  public boolean integrity() {
    for (int i = 0; i < mScores.length; i++) {
      Exam.assertEquals(mScores[i].length, mScores.length);
      for (int j = 0; j < mScores[i].length; j++) {
        final int sc = mScores[i][j];
        Exam.assertTrue(-10 < sc && sc < 20);
        Exam.assertEquals(sc, mScores[j][i]);
      }
    }
    return true;
  }

  @Override
  public boolean globalIntegrity() {
    return integrity();
  }

  protected int codeIndex(final ProteinFastaSymbolTable table, final String x) {
    //final int y;
    final Protein pr = (Protein) table.scanResidue(x.charAt(0));
//    if (pr.ignore() && !x.equals("X") && !x.equals("*")) {
//      y = -1;
//    } else {
//      y = pr.ordinal();
//    }
    return pr.ordinal();
  }

  protected void parse(final BufferedReader in) throws IOException {
    final ProteinFastaSymbolTable table = new ProteinFastaSymbolTable();
    int[] colids = null;
    String line;
    while ((line = in.readLine()) != null) {
      if (!line.startsWith("#")) {
        //get the ids for the columns
        final String[] splits = line.split("\\ +");
        colids = new int[splits.length];
        for (int i = 1; i < splits.length; i++) {
          final String x = splits[i];
          colids[i] = codeIndex(table, x);
        }
        break;
      }
    }
    //System.err.println(Arrays.toString(colids));
    while ((line = in.readLine()) != null) {
      final String[] splits = line.split("\\ +");
      final int row = codeIndex(table, splits[0]);
      if (row > -1) {
        //System.err.println("row=" + row);
        for (int i = 1; i < splits.length; i++) {
          final int col = colids[i];
          if (col == -1) {
            continue;
          }
          final int sc = Integer.parseInt(splits[i]);
          mScores[row][col] = sc;
          //System.err.println("[" + i + "]" + splits[i].length() + ":" + splits[i]);
        }
      }
    }
  }
}
