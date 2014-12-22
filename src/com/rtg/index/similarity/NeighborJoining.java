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

import java.util.ArrayList;

import com.rtg.util.PortableRandom;
import com.rtg.util.diagnostic.Diagnostic;


/**
 * An implementation of the Neighbor-joining algorithm of Saitou and Nei.
 * This implementation is based on the description in Wikipedia.
 *
 */
public class NeighborJoining {

  private final PortableRandom mRandom;

  /**
   * Creates a new Neighbor Joining object.
   */
  public NeighborJoining() {
    this(23);
  }

  /**
   * Creates a new Neighbor Joining object.
   * @param seed random number generator seed
   */
  public NeighborJoining(final long seed) {
    mRandom = new PortableRandom(seed);
  }

  /* Retrieve the distance between two taxa. */
  private static double distance(final int i, final int j, final ArrayList<ArrayList<Double>> d) {
    return i == j ? 0 : (i > j ? d.get(i).get(j) : d.get(j).get(i));
  }

  /* Remove a given item from the distance matrix. */
  private static void remove(final int f, final ArrayList<BinaryTree> s, final ArrayList<ArrayList<Double>> d) {
    for (int k = f + 1; k < d.size(); k++) {
      d.get(k).remove(f);
    }
    s.remove(f);
    d.remove(f);
  }

  /**
   * Perform the neighbor-joining algorithm.  Given a list of node names and a
   * matrix of distances between sequences, infer a binary tree
   * using the neighbor-joining algorithm.
   *
   * @param nodeNames the list of node names for translating back from the matrix.
   * @param matrix count of number of hashes in common between pairs of sequences.
   * @return the resulting tree
   */
  public BinaryTree neighborJoin(final ArrayList<String> nodeNames, final SimilarityMatrix matrix) {
    final ArrayList<ArrayList<Double>> d = makeArray(matrix);
    return neighborJoin(nodeNames, d);
  }

  /**
   * Produce an array which is a lower triangular matrix with normalized and inverted
   * similarity scores (the matrix contains counts of common hashes for pairs of sequences).
   * @param matrix with the original counts.
   * @return the lower triangular array.
   */
  static ArrayList<ArrayList<Double>> makeArray(final SimilarityMatrix matrix) {
    final ArrayList<ArrayList<Double>> d;
    final int length = matrix.length();
    //compute normalization factors
    final double[] norm = new double[length];
    for (int i = 0; i < length; i++) {
      //be careful to deal with 0's on diagonal
      final double n = matrix.get(i, i);
      assert n >= 0;
      norm[i] = Math.sqrt(n == 0 ? 1 : n);
    }
    d = new ArrayList<>();
    for (int i = 0; i < length; i++) {
      final ArrayList<Double> row = new ArrayList<>();
      for (int j = 0; j < i; j++) {
        final double v = matrix.get(i, j) / (norm[i] * norm[j]);
        assert v >= 0.0 && !Double.isInfinite(v) && !Double.isNaN(v) : v;
        final double w = 1.0 / (1.0 + v);
        assert w >= 0.0 && !Double.isInfinite(w) && !Double.isNaN(w) : w;
        row.add(w);
      }
      d.add(row);
    }
    return d;
  }

  /**
   * Perform the neighbor-joining algorithm.  Given a list of node names and a
   * lower triangular matrix of distances between names, infer a binary tree
   * using the neighbor-joining algorithm.
   *
   * @param nodeNames the node names
   * @param d distance matrices
   * @return the resulting tree
   */
  BinaryTree neighborJoin(final ArrayList<String> nodeNames, final ArrayList<ArrayList<Double>> d) {
    final ArrayList<BinaryTree> s = new ArrayList<>();
    for (final String name : nodeNames) {
      s.add(new BinaryTree(null, null, 0, 0, name));
    }
    while (s.size() > 1) {
      assert d.size() == s.size();
      // Precompute column sums
      final double[] colSum = new double[d.size()];
      for (int j = 0; j < d.size(); j++) {
        double sum = 0;
        for (int k = 0; k < d.size(); k++) {
          sum += distance(j, k, d);
        }
        colSum[j] = sum;
      }

      // Compute minimal entry in (implicit) Q matrix
      final int r = d.size() - 2;
      double best = Double.POSITIVE_INFINITY;
      int c = 0;
      int f = -1;
      int g = -1;
      for (int k = 1; k < d.size(); k++) {
        for (int j = 0; j < k; j++) {
          final double q = r * distance(k, j, d) - colSum[k] - colSum[j];
          if (Double.doubleToRawLongBits(q) == Double.doubleToRawLongBits(best) && mRandom.nextInt(++c) == 0) {
            // Break ties fairly
            f = k;
            g = j;
          } else if (q < best) {
            c = 1;
            f = k;
            g = j;
            best = q;
          }
        }
      }
      assert f >= 0 && g >= 0 : "f=" + f + " g=" + g;

      // Compute distance of merged node to new node
      final double dfu, dgu, dfg = 0.5 * distance(f, g, d);
      if (d.size() > 2) {
        dfu = dfg + 0.5 * (colSum[f] - colSum[g]) / (d.size() - 2);
        dgu = dfg + 0.5 * (colSum[g] - colSum[f]) / (d.size() - 2);
      } else {
        assert Math.abs(colSum[f] - colSum[g]) < 0.0000000001;
        dfu = dfg;
        dgu = dfg;
      }

      // Compute distance of merged node to all other nodes
      final ArrayList<Double> newRow = new ArrayList<>();
      for (int i = 0; i < d.size(); i++) {
        if (i != f && i != g) {
          newRow.add(0.5 * (distance(f, i, d) - dfu + distance(g, i, d) - dgu));
        }
      }
      final BinaryTree newNode = new BinaryTree(s.get(f), s.get(g), dfu, dgu, String.valueOf(d.size()));
      Diagnostic.userLog("NeighborJoining: " + s.get(f).getLabel() + " + " + s.get(g).getLabel() + " -> " + newNode.getLabel());

      // Remove merged nodes starting with the larger one
      if (f > g) {
        remove(f, s, d);
        remove(g, s, d);
      } else {
        remove(g, s, d);
        remove(f, s, d);
      }

      // Add combined node at end of name list and matrix
      assert newRow.size() == d.size();
      d.add(newRow);
      s.add(newNode);
    }
    return s.get(0);
  }
}

