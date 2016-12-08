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
package com.rtg.index.similarity;

import java.io.IOException;
import java.util.List;

import com.rtg.util.StringUtils;
import com.rtg.util.Utils;
import com.rtg.util.diagnostic.ErrorType;
import com.rtg.util.diagnostic.SlimException;
import com.rtg.util.integrity.Exam;
import com.rtg.util.integrity.IntegralAbstract;

/**
 * Accumulates values which indicate the similarity between two sequences.
 * This version is used only for the case of all sequences compared against themselves.
 */
public class SimilarityMatrix extends IntegralAbstract {

  private final int mLength;

  private final double[][] mCounts;


  /**
   * @param numberSequences same as size of matrix.
   */
  public SimilarityMatrix(final long numberSequences) {
    if (numberSequences > Integer.MAX_VALUE) {
      throw new RuntimeException(numberSequences + "");
    }
    mLength = (int) numberSequences;
    mCounts = new double[mLength][];
    for (int i = 0; i < mLength; ++i) {
      mCounts[i] = new double[i + 1];
    }
  }

  /**
   * Get the size of the array.
   * @return the size of the array (it is a length x length array).
   */
  public int length() {
    return mLength;
  }

  /**
   * Increment count taking into account commutativity of matrix.
   * @param a first index.
   * @param b second index.
   */
  public void increment(final int a, final int b) {
    if (a < b) {
      mCounts[b][a]++;
    } else {
      mCounts[a][b]++;
    }
  }

  /**
   * Increment count taking into account commutativity of matrix.
   * @param a first index.
   * @param b second index.
   * @param incr amount to increment count by.
   */
  public void increment(final int a, final int b, final double incr) {
    assert incr > 0;
    if (a < b) {
      mCounts[b][a] += incr;
    } else {
      mCounts[a][b] += incr;
    }
  }

  /**
   * Get count taking into account commutativity of matrix.
   * @param a first index.
   * @param b second index.
   * @return the count (&gt;= 0).
   */
  public double get(final int a, final int b) {
    if (a < b) {
      return mCounts[b][a];
    } else {
      return mCounts[a][b];
    }
  }

  /**
   * Set count taking into account commutativity of matrix.
   * Provided for testing when need mock matrices.
   * @param a first index.
   * @param b second index.
   * @param v the value to be assigned.
   */
  public void set(final int a, final int b, final double v) {
    if (v < 0) {
      throw new RuntimeException();
    }
    if (a < b) {
      mCounts[b][a] = v;
    } else {
      mCounts[a][b] = v;
    }
  }

  @Override
  public void toString(final StringBuilder sb) {
    sb.append("SimilarityMatrix ").append(mLength).append(StringUtils.LS);
    for (int i = 0; i < mLength; ++i) {
      sb.append("[").append(i).append("]");
      for (int j = 0; j < mLength; ++j) {
        sb.append("\t");
        sb.append(Utils.realFormat(get(i, j), 0));
      }
      sb.append(StringUtils.LS);
    }
    sb.append(StringUtils.LS);
  }

  /**
   * Write the contents of this matrix to the given appendable.
   * @param out place to write to.
   * @param names the names of the sequences.
   * @throws IOException if an I/O error occurs.
   */
  public void write(final Appendable out, List<String> names) throws IOException {
    if (names.size() != mLength) {
      throw new SlimException(ErrorType.INFO_ERROR, "Size of matrix not equal to number of labels provided.");
    }
    final String tab = "\t";
    for (int i = 0; i < mLength; ++i) {
      out.append(tab);
      out.append(names.get(i));
    }
    out.append(StringUtils.LS);
    for (int i = 0; i < mLength; ++i) {
      out.append(names.get(i));
      for (int j = 0; j < mLength; ++j) {
        out.append(tab);
        out.append(Utils.realFormat(get(i, j), 0));
      }
      out.append(StringUtils.LS);
    }
  }

  /**
   * Function to turn the matrix into a printable string with appropriate name labels.
   * @param names the names of the sequences.
   * @return the printable matrix.
   */
  public String toString(final List<String> names) {
    final StringBuilder sb = new StringBuilder();
    try {
      write(sb, names);
    } catch (final IOException e) {
      throw new RuntimeException(e);
    }
    return sb.toString();
  }

  @Override
  public boolean globalIntegrity() {
    integrity();
    for (int i = 0; i < mLength; ++i) {
      Exam.assertEquals(i + 1, mCounts[i].length);
      for (int j = 0; j <= i; ++j) {
        Exam.assertTrue(mCounts[i][j] >= 0);
      }
    }
    return true;
  }

  @Override
  public boolean integrity() {
    if (mCounts == null) {
      throw new NullPointerException();
    } else {
      for (int i = 0; i < mLength; ++i) {
        Exam.assertEquals(i + 1, mCounts[i].length);
      }
    }
    return true;
  }
}
