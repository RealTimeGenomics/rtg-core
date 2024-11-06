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
package com.rtg.metagenomics.matrix;

/**
 */
public final class MatrixUtils {

  private MatrixUtils() { } //prevent instantiation

  /**
   * Return the transpose of a matrix.
   * @param a matrix to be transposed.
   * @return the transposed version (may be the original version a).
   */
  public static Matrix transpose(final Matrix a) {
    if (a.isSymmetric()) {
      return a;
    }
    return new MatrixTranspose(a);
  }

  /**
   * Compute A B.
   * @param a matrix A.
   * @param b matrix B.
   * @return result of multiplication.
   */
  public static Matrix multiply(final Matrix a, final Matrix b) {
    if (a.cols() != b.rows()) {
      throw new IllegalArgumentException("A.length=" + a.cols() + " B.length=" + b.rows());
    }
    final Matrix res = new MatrixSimple(a.rows(), b.cols());
    for (int i = 0; i < a.rows(); ++i) {
      for (int j = 0; j < b.cols(); ++j) {
        double sum = 0.0;
        for (int k = 0; k < a.cols(); ++k) {
          sum += a.get(i, k) * b.get(k, j);
        }
        res.set(i, j, sum);
      }
    }
    return res;
  }

  /**
   * Compute A v.
   * @param a matrix A.
   * @param v vector.
   * @return result of multiplication.
   */
  public static Vector multiply(final Matrix a, final Vector v) {
    final int n = a.rows();
    if (n != v.size()) {
      throw new IllegalArgumentException("A.length=" + n + " v.length=" + v.size());
    }
    final Vector res = new Vector(n);
    for (int i = 0; i < n; ++i) {
      double sum = 0.0;
      for (int j = 0; j < n; ++j) {
        sum += v.get(j) * a.get(i, j);
      }
      res.set(i, sum);
    }
    return res;
  }

  /**
   * Multiply a vector by a constant.
   * @param a constant.
   * @param v vector to be multiplied.
   * @return vector result.
   */
  public static Vector multiply(final double a, final Vector v) {
    final int n = v.size();
    final Vector res = new Vector(n);
    for (int i = 0; i < n; ++i) {
      res.set(i, a * v.get(i));
    }
    return res;
  }

  /**
   * Compute dot product of two vectors.
   * @param v first vector.
   * @param w second vector.
   * @return dot product.
   */
  public static double multiply(final Vector v, final Vector w) {
    final int n = v.size();
    if (n != w.size()) {
      throw new IllegalArgumentException("v.length=" + n + " w.length=" + w.size());
    }
    double sum = 0.0;
    for (int i = 0; i < n; ++i) {
      sum += v.get(i) * w.get(i);
    }
    return sum;
  }

  /**
   * Compute the square of the norm of a vector.
   * @param v vector whose norm is to be computed.
   * @return  the square of the norm.
   */
  public static double norm2(final Vector v) {
    return multiply(v, v);
  }

  /**
   * Add two vectors.
   * @param v first vector.
   * @param w second vector.
   * @return the vector sum.
   */
  public static Vector add(final Vector v, final Vector w) {
    final int n = v.size();
    if (n != w.size()) {
      throw new IllegalArgumentException("v.length=" + n + " w.length=" + w.size());
    }
    final Vector res = new Vector(n);
    for (int i = 0; i < n; ++i) {
      res.set(i, v.get(i) + w.get(i));
    }
    return res;
  }

  /**
   * Subtract two vectors.
   * @param v first vector.
   * @param w second vector.
   * @return the vector difference (v - w).
   */
  public static Vector subtract(final Vector v, final Vector w) {
    final int n = v.size();
    if (n != w.size()) {
      throw new IllegalArgumentException("v.length=" + n + " w.length=" + w.size());
    }
    final Vector res = new Vector(n);
    for (int i = 0; i < n; ++i) {
      res.set(i, v.get(i) - w.get(i));
    }
    return res;
  }

  /**
   * Compute the inner product, <code>v A w</code>, generated by the matrix A.
   * <code>r = sum_{i,j} v_i A_{i,j} w_j</code>
   * @param v first vector.
   * @param a matrix A.
   * @param w second vector.
   * @return the inner product.
   */
  public static double innerProduct(final Vector v, final Matrix a, final Vector w) {
    double sum = 0.0;
    final int n = a.rows();
    if (n != v.size() || n != w.size()) {
      throw new IllegalArgumentException("v.length=" + v.size() + " A.length=" + a.rows() + " w.length=" + w.size());
    }
    for (int i = 0; i < n; ++i) {
      final double vi = v.get(i);
      for (int j = 0; j < n; ++j) {
        sum += vi * a.get(i, j) * w.get(j);
      }
    }
    return sum;
  }

  /**
   * Construct a new matrix with the same length and symmetry properties as the original.
   * @param a matrix to be copied (but not the values).
   * @return a new matrix with the same shape as <code>a</code>.
   */
  public static Matrix sameShape(final Matrix a) {
    final int n = a.rows();
    final Matrix res;
    if (a.isSymmetric()) {
      res = new MatrixSymmetric(n);
    } else {
      res = new MatrixSimple(n);
    }
    return res;
  }

  /**
   * Compute the point product, <code>v A w</code>, generated by the matrix A.
   * <code>R_{i,j} = v_i A_{i,j} w_j</code>
   *
   * @param v first vector.
   * @param a matrix A.
   * @param w second vector.
   * @return the inner product.
   */
  public static Matrix pointProduct(final Vector v, final Matrix a, final Vector w) {
    final Matrix res = new MatrixSimple(a.rows());
    final int n = a.rows();
    if (n != v.size() || n != w.size()) {
      throw new IllegalArgumentException("v.length=" + v.size() + " A.length=" + a.rows() + " w.length=" + w.size());
    }
    for (int i = 0; i < n; ++i) {
      final double vi = v.get(i);
      for (int j = 0; j < n; ++j) {
        final double s = vi * a.get(i, j) * w.get(j);
        res.set(i, j, s);
      }
    }
    return res;
  }

  /**
   * Compute the point product, <code>v A v^T</code>, generated by the matrix A.
   * <code>R_{i,j} = v_i A_{i,j} v_j</code>
   *
   * @param v first vector.
   * @param a matrix A.
   * @return the inner product.
   */
  public static Matrix pointProduct(final Vector v, final Matrix a) {
    //if a is symmetrix then so is the result and the calculation can be optimized
    final Matrix res = sameShape(a);
    final int n = a.rows();
    if (n != v.size()) {
      throw new IllegalArgumentException("v.length=" + v.size() + " A.length=" + a.rows());
    }
    for (int i = 0; i < n; ++i) {
      final double vi = v.get(i);
      final int hi = a.isSymmetric() ? i + 1 : n;
      for (int j = 0; j < hi; ++j) {
        final double s = vi * a.get(i, j) * v.get(j);
        res.set(i, j, s);
      }
    }
    return res;
  }

  /**
   * Compute the point product, <code>v w</code>.
   * <code>R_i = v_i  w_i</code>
   *
   * @param v first vector.
   * @param w second vector.
   * @return the inner product.
   */
  public static Vector pointProduct(final Vector v, final Vector w) {
    final int n = v.size();
    if (n != w.size()) {
      throw new IllegalArgumentException("v.length=" + v.size() + " w.length=" + w.size());
    }
    final Vector res = new Vector(n);
    for (int i = 0; i < n; ++i) {
      res.set(i, v.get(i) * w.get(i));
    }
    return res;
  }

  /**
   * Extract the trace (diagonal elements from a matrix).
   * @param a matrix from which the trace is to be extracted.
   * @return a vector containing the trace.
   */
  public static Vector trace(final Matrix a) {
    final int n = a.rows();
    final Vector res = new Vector(n);
    for (int i = 0; i < n; ++i) {
      res.set(i, a.get(i, i));
    }
    return res;
  }

  /**
   * Compute the point by point inverse of a vector.
   * @param v vector to be inverted.
   * @return a vector with the inverse.
   */
  public static Vector inverse(final Vector v) {
    final int n = v.size();
    final Vector res = new Vector(n);
    for (int i = 0; i < n; ++i) {
      res.set(i, 1.0 / v.get(i));
    }
    return res;
  }

  /**
   * Compute the point by point square root of a vector.
   * @param v vector to be used.
   * @return a vector with the square roots.
   */
  public static Vector sqrt(final Vector v) {
    final int n = v.size();
    final Vector res = new Vector(n);
    for (int i = 0; i < n; ++i) {
      res.set(i, Math.sqrt(v.get(i)));
    }
    return res;
  }

  /**
   * Compute the negation of a vector.
   * @param v vector to be negated.
   * @return negated vector.
   */
  public static Vector negative(final Vector v) {
    final int dim = v.size();
    final Vector neg = new Vector(dim);
    for (int i = 0; i < dim; ++i) {
      final double d = -v.get(i);
      neg.set(i, d);
    }
    return neg;
  }

  /**
   * Check that all the values in a vector are finite.
   * @param v vector to be checked.
   * @return true iff all values in the vector are finite.
   */
  public static boolean isFinite(final Vector v) {
    for (int i = 0; i < v.size(); ++i) {
      final double vv = v.get(i);
      if (Double.isNaN(vv) || Double.isInfinite(vv)) {
        return false;
      }
    }
    return true;
  }

  /**
   * Check that all the values in a vector are finite and positive.
   * @param v vector to be checked.
   * @return true iff all values in the vector are finite and positive.
   */
  public static boolean isFinitePositive(final Vector v) {
    for (int i = 0; i < v.size(); ++i) {
      final double vv = v.get(i);
      if (vv < 0.0 || Double.isNaN(vv) || Double.isInfinite(vv)) {
        return false;
      }
    }
    return true;
  }

  /**
   * Sum all values in a vector.
   * @param v vector to be summed.
   * @return trace
   */
  public static double trace(final Vector v) {
    double sum = 0.0;
    for (int i = 0; i < v.size(); ++i) {
      final double vv = v.get(i);
      sum += vv;
    }
    return sum;
  }
}
