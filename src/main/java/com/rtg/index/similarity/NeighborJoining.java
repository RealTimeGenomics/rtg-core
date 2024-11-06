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
package com.rtg.index.similarity;

import java.util.ArrayList;
import java.util.List;

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
  private static double distance(final int i, final int j, final List<List<Double>> d) {
    return i == j ? 0 : (i > j ? d.get(i).get(j) : d.get(j).get(i));
  }

  /* Remove a given item from the distance matrix. */
  private static void remove(final int f, final List<BinaryTree> s, final List<List<Double>> d) {
    for (int k = f + 1; k < d.size(); ++k) {
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
  public BinaryTree neighborJoin(final List<String> nodeNames, final SimilarityMatrix matrix) {
    final List<List<Double>> d = makeArray(matrix);
    return neighborJoin(nodeNames, d);
  }

  /**
   * Produce an array which is a lower triangular matrix with normalized and inverted
   * similarity scores (the matrix contains counts of common hashes for pairs of sequences).
   * @param matrix with the original counts.
   * @return the lower triangular array.
   */
  static List<List<Double>> makeArray(final SimilarityMatrix matrix) {
    final List<List<Double>> d;
    final int length = matrix.length();
    //compute normalization factors
    final double[] norm = new double[length];
    for (int i = 0; i < length; ++i) {
      //be careful to deal with 0's on diagonal
      final double n = matrix.get(i, i);
      assert n >= 0;
      norm[i] = Math.sqrt(n == 0 ? 1 : n);
    }
    d = new ArrayList<>();
    for (int i = 0; i < length; ++i) {
      final List<Double> row = new ArrayList<>();
      for (int j = 0; j < i; ++j) {
        final double v = matrix.get(i, j) / (norm[i] * norm[j]);
        assert v >= 0.0 && Double.isFinite(v) : v;
        final double w = 1.0 / (1.0 + v);
        assert w >= 0.0 && Double.isFinite(w) : w;
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
  BinaryTree neighborJoin(final List<String> nodeNames, final List<List<Double>> d) {
    final List<BinaryTree> s = new ArrayList<>();
    for (final String name : nodeNames) {
      s.add(new BinaryTree(null, null, 0, 0, name));
    }
    while (s.size() > 1) {
      assert d.size() == s.size() : "distance size=" + d.size() + " cf. names size=" + s.size();
      // Precompute column sums
      final double[] colSum = new double[d.size()];
      for (int j = 0; j < d.size(); ++j) {
        double sum = 0;
        for (int k = 0; k < d.size(); ++k) {
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
      for (int k = 1; k < d.size(); ++k) {
        for (int j = 0; j < k; ++j) {
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
      for (int i = 0; i < d.size(); ++i) {
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

