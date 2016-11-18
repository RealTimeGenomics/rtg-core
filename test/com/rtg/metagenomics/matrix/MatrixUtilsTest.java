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
package com.rtg.metagenomics.matrix;

import static com.rtg.util.StringUtils.LS;

import com.rtg.util.integrity.Exam;

import junit.framework.TestCase;

/**
 */
public class MatrixUtilsTest extends TestCase {

  public void testTranspose1() {
    final Matrix a = new MatrixSymmetric(2);
    a.set(0, 0, 4.0);
    a.set(1, 1, 3.0);
    a.set(0, 1, 1.0);
    final Matrix mt = MatrixUtils.transpose(a);
    assertEquals(4.0, mt.get(0, 0));
    assertEquals(3.0, mt.get(1, 1));
    assertEquals(1.0, mt.get(0, 1));
    assertEquals(1.0, mt.get(1, 0));
  }

  public void testTranspose2() {
    final Matrix a = new MatrixSimple(2);
    a.set(0, 0, 4.0);
    a.set(1, 1, 3.0);
    a.set(0, 1, 1.0);
    final Matrix mt = MatrixUtils.transpose(a);
    assertEquals(4.0, mt.get(0, 0));
    assertEquals(3.0, mt.get(1, 1));
    assertEquals(0.0, mt.get(0, 1));
    assertEquals(1.0, mt.get(1, 0));
  }

  public void testMultiplyMatrixMatrix() {
    final Matrix a = new MatrixSimple(new double[][] {{0.5, 1.0}, {0.2, 1.0}});
    final Matrix b = new MatrixSimple(new double[][] {{1.0, 3.0}, {0.0, 2.0}});
    final Matrix c = MatrixUtils.multiply(a, b);
    final String str = c.toString();
    //System.err.println(str);
    final String exp = ""
      + "[0]  0.5000  3.5000" + LS
      + "[1]  0.2000  2.6000" + LS
      ;
    assertEquals(exp, str);
  }

  public final void testMultiplyMatrixBad() {
    final Matrix a = new MatrixSymmetric(2);
    final Matrix v = new MatrixSymmetric(1);
    try {
      MatrixUtils.multiply(a, v);
      fail();
    } catch (final IllegalArgumentException e) {
      //expected
      assertEquals("A.length=2 B.length=1", e.getMessage());
    }
  }


  public final void testMultiplyMatrixVector() {
    final Matrix a = new MatrixSymmetric(2);
    a.set(0, 0, 4.0);
    a.set(1, 1, 3.0);
    a.set(0, 1, 1.0);
    final Vector v = new Vector(new double[] {2.0, 1.0});
    final Vector w = MatrixUtils.multiply(a, v);
    Exam.assertEquals(2, w.size());
    Exam.assertEquals(9.0,  w.get(0));
    Exam.assertEquals(5.0, w.get(1));
  }

  public final void testMultiplyMatrixVectorBad() {
    final Matrix a = new MatrixSymmetric(2);
    final Vector v = new Vector(1);
    try {
      MatrixUtils.multiply(a, v);
      fail();
    } catch (final IllegalArgumentException e) {
      //expected
      assertEquals("A.length=2 v.length=1", e.getMessage());
    }
  }

  public final void testMultiplyDoubleVector0() {
    final Vector v = new Vector(0);
    final Vector w = MatrixUtils.multiply(2.0, v);
    Exam.assertEquals(0, w.size());
  }

  public final void testMultiplyDoubleVector() {
    final Vector v = new Vector(new double[] {0.3, -1.5});
    final Vector w = MatrixUtils.multiply(2.0, v);
    Exam.assertEquals(2, w.size());
    Exam.assertEquals(0.6,  w.get(0));
    Exam.assertEquals(-3.0, w.get(1));
  }

  public final void testMultiplyVectorVector() {
    final Vector v = new Vector(new double[] {0.3, -1.5});
    final Vector w = new Vector(new double[] {-0.3, -1.0});
    Exam.assertEquals(1.41, MatrixUtils.multiply(v, w));
  }

  public final void testMultiplyVectorVectorBad() {
    final Vector v = new Vector(2);
    final Vector w = new Vector(1);
    try {
      MatrixUtils.multiply(v, w);
      fail();
    } catch (final IllegalArgumentException e) {
      //expected
      assertEquals("v.length=2 w.length=1", e.getMessage());
    }
  }

  public final void testNorm2() {
    final Vector v = new Vector(new double[] {0.3, -1.5});
    Exam.assertEquals(2.34, MatrixUtils.norm2(v));
  }

  public final void testAdd() {
    final Vector v = new Vector(new double[] {0.3, -1.5});
    final Vector w = new Vector(new double[] {0.8, 0.5});
    final Vector x = MatrixUtils.add(v, w);
    Exam.assertEquals(2, x.size());
    Exam.assertEquals(1.1,  x.get(0));
    Exam.assertEquals(-1.0, x.get(1));
  }

  public final void testAddBad() {
    final Vector v = new Vector(2);
    final Vector w = new Vector(1);
    try {
      MatrixUtils.add(v, w);
      fail();
    } catch (final IllegalArgumentException e) {
      //expected
      assertEquals("v.length=2 w.length=1", e.getMessage());
    }
  }

  public final void testSubtract() {
    final Vector v = new Vector(new double[] {0.3, -1.5});
    final Vector w = new Vector(new double[] {0.8, 0.5});
    final Vector x = MatrixUtils.subtract(v, w);
    Exam.assertEquals(2, x.size());
    Exam.assertEquals(-0.5,  x.get(0));
    Exam.assertEquals(-2.0, x.get(1));
  }

  public final void testSubtractBad() {
    final Vector v = new Vector(2);
    final Vector w = new Vector(1);
    try {
      MatrixUtils.subtract(v, w);
      fail();
    } catch (final IllegalArgumentException e) {
      //expected
      assertEquals("v.length=2 w.length=1", e.getMessage());
    }
  }

  public final void testInnerProduct() {
    final Matrix a = new MatrixSymmetric(2);
    a.set(0, 0, 4.0);
    a.set(1, 1, 3.0);
    a.set(0, 1, 1.0);
    final Vector v = new Vector(new double[] {2.0, 1.0});
    final Vector w = new Vector(new double[] {0.5, 1.0});
    final double x = MatrixUtils.innerProduct(v, a, w);
    Exam.assertEquals(9.5, x);
  }

  public final void testInnerProductBad() {
    final Matrix a = new MatrixSymmetric(2);
    final Vector v = new Vector(2);
    final Vector w = new Vector(1);
    try {
      MatrixUtils.innerProduct(v, a, w);
      fail();
    } catch (final IllegalArgumentException e) {
      //expected
      assertEquals("v.length=2 A.length=2 w.length=1", e.getMessage());
    }
    try {
      MatrixUtils.innerProduct(w, a, v);
      fail();
    } catch (final IllegalArgumentException e) {
      //expected
      assertEquals("v.length=1 A.length=2 w.length=2", e.getMessage());
    }
  }

  public final void testSameShapeSimple() {
    final Matrix a = new MatrixSimple(3);
    assertEquals(3, a.size());
    assertTrue(a.isSymmetric());
    a.set(1, 2, Math.PI);
    assertFalse(a.isSymmetric());
    final Matrix b = MatrixUtils.sameShape(a);
    assertEquals(3, b.size());
    assertTrue(b.isSymmetric());
  }

  public final void testSameShapeSym() {
    final Matrix a = new MatrixSymmetric(3);
    assertEquals(3, a.size());
    assertTrue(a.isSymmetric());
    final Matrix b = MatrixUtils.sameShape(a);
    assertEquals(3, b.size());
    assertTrue(b.isSymmetric());
  }

  public final void testPointProduct() {
    final Matrix a = new MatrixSymmetric(2);
    a.set(0, 0, 4.0);
    a.set(1, 1, 3.0);
    a.set(0, 1, 1.0);
    final Vector v = new Vector(new double[] {2.0, 1.0});
    final Vector w = new Vector(new double[] {0.33, 1.0});
    final Matrix x = MatrixUtils.pointProduct(v, a, w);
    assertFalse(x.isSymmetric());
    final String exp = ""
      + "[0]  2.6400  2.0000" + LS
      + "[1]  0.3300  3.0000" + LS
      ;
    assertEquals(exp, x.toString());
  }

  public final void testPointProductBad() {
    final Matrix a = new MatrixSymmetric(2);
    final Vector v = new Vector(2);
    final Vector w = new Vector(1);
    try {
      MatrixUtils.pointProduct(v, a, w);
      fail();
    } catch (final IllegalArgumentException e) {
      //expected
      assertEquals("v.length=2 A.length=2 w.length=1", e.getMessage());
    }
    try {
      MatrixUtils.pointProduct(w, a, v);
      fail();
    } catch (final IllegalArgumentException e) {
      //expected
      assertEquals("v.length=1 A.length=2 w.length=2", e.getMessage());
    }
  }

  public final void testPointProductVV() {
    final Vector v = new Vector(new double[] {2.0, 1.0});
    final Vector w = new Vector(new double[] {0.33, 1.0});
    final Vector x = MatrixUtils.pointProduct(v, w);
    assertEquals(2, x.size());
    assertEquals(0.66, x.get(0));
    assertEquals(1.00, x.get(1));

  }

  public final void testPointProductVVBad() {
    final Vector v = new Vector(2);
    final Vector w = new Vector(1);
    try {
      MatrixUtils.pointProduct(v, w);
      fail();
    } catch (final IllegalArgumentException e) {
      //expected
      assertEquals("v.length=2 w.length=1", e.getMessage());
    }
  }

  public final void testPointProduct1() {
    final Matrix a = new MatrixSymmetric(2);
    a.set(0, 0, 4.0);
    a.set(1, 1, 3.0);
    a.set(0, 1, 1.0);
    final Vector v = new Vector(new double[] {2.0, 1.0});
    final Matrix x = MatrixUtils.pointProduct(v, a);
    assertTrue(x.isSymmetric());
    final String exp = ""
      + "[0]  16.0000" + LS
      + "[1]  2.0000  3.0000" + LS
      ;
    assertEquals(exp, x.toString());
  }

  public final void testPointProduct2() {
    final Matrix a = new MatrixSimple(new double[][] {{4.0, 3.0}, {3.0, 1.0}});
    final Vector v = new Vector(new double[] {2.0, 1.0});
    final Matrix x = MatrixUtils.pointProduct(v, a);
    assertTrue(x.isSymmetric());
    final String exp = ""
      + "[0]  16.0000" + LS
      + "[1]  6.0000  1.0000" + LS
      ;
    assertEquals(exp, x.toString());
  }

  public final void testPointProductBad1() {
    final Matrix a = new MatrixSymmetric(2);
    final Vector w = new Vector(1);
    try {
      MatrixUtils.pointProduct(w, a);
      fail();
    } catch (final IllegalArgumentException e) {
      //expected
      assertEquals("v.length=1 A.length=2", e.getMessage());
    }
  }

  public void testTrace() {
    final Matrix a = new MatrixSymmetric(2);
    a.set(0, 0, 4.0);
    a.set(1, 1, 16.0);
    a.set(0, 1, 1.0);
    final Vector v = MatrixUtils.trace(a);
    assertEquals(2, v.size());
    assertEquals(4.0, v.get(0));
    assertEquals(16.0, v.get(1));
  }

  public void testInverse() {
    final Vector v = new Vector(new double[]{1.0, 2.0});
    final Vector w = MatrixUtils.inverse(v);
    assertEquals(2, w.size());
    assertEquals(1.0, w.get(0));
    assertEquals(0.5, w.get(1));
  }

  public void testSqrt() {
    final Vector v = new Vector(new double[]{4.0, 9.0});
    final Vector w = MatrixUtils.sqrt(v);
    assertEquals(2, w.size());
    assertEquals(2.0, w.get(0));
    assertEquals(3.0, w.get(1));
  }

  public void testNegative() {
    final Vector v = new Vector(new double[]{4.0, -9.0});
    final Vector w = MatrixUtils.negative(v);
    assertEquals(2, w.size());
    assertEquals(-4.0, w.get(0));
    assertEquals(9.0, w.get(1));
  }

  public void testIsFinite() {
    assertTrue(MatrixUtils.isFinite(new Vector(new double[]{4.0, -9.0})));
    assertTrue(MatrixUtils.isFinite(new Vector(new double[]{})));

    assertFalse(MatrixUtils.isFinite(new Vector(new double[]{Double.NaN, -9.0})));
    assertFalse(MatrixUtils.isFinite(new Vector(new double[]{4.0, Double.POSITIVE_INFINITY})));
    assertFalse(MatrixUtils.isFinite(new Vector(new double[]{4.0, Double.NEGATIVE_INFINITY})));
  }

  public void testIsFinitePositive() {
    assertTrue(MatrixUtils.isFinitePositive(new Vector(new double[]{4.0, 0.0})));
    assertTrue(MatrixUtils.isFinitePositive(new Vector(new double[]{})));

    assertFalse(MatrixUtils.isFinitePositive(new Vector(new double[]{4.0, -9.0})));
    assertFalse(MatrixUtils.isFinitePositive(new Vector(new double[]{Double.NaN, 4.0})));
    assertFalse(MatrixUtils.isFinitePositive(new Vector(new double[]{4.0, Double.POSITIVE_INFINITY})));
    assertFalse(MatrixUtils.isFinitePositive(new Vector(new double[]{4.0, Double.NEGATIVE_INFINITY})));
  }
}
