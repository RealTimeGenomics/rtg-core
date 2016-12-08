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
package com.rtg.util.memory;

import java.lang.reflect.Field;
import java.util.Iterator;

import junit.framework.TestCase;

/**
 */
public class ClassMemoryTest extends TestCase {

  public void testBooleanPrim() {
    primTst(Boolean.TYPE, MemoryUsage.BOOLEAN_SIZE);
  }

  public void testBytePrim() {
    primTst(Byte.TYPE, MemoryUsage.BYTE_SIZE);
  }

  public void testCharPrim() {
    primTst(Character.TYPE, MemoryUsage.CHAR_SIZE);
  }

  public void testShortPrim() {
    primTst(Short.TYPE, MemoryUsage.SHORT_SIZE);
  }

  public void testIntPrim() {
    primTst(Integer.TYPE, MemoryUsage.INT_SIZE);
  }

  public void testFloatPrim() {
    primTst(Float.TYPE, MemoryUsage.FLOAT_SIZE);
  }

  public void testLongPrim() {
    primTst(Long.TYPE, MemoryUsage.LONG_SIZE);
  }

  public void testDoublePrim() {
    primTst(Double.TYPE, MemoryUsage.DOUBLE_SIZE);
  }

  private static void primTst(final Class<?> classId, final int size) {
    final ClassMemory.Info info = ClassMemory.getMemoryInfo(classId);
    assertTrue(info.isPrimitive());

    assertEquals(info + "", size, info.getSize());

    final int arr = 0;
    assertEquals(arr, info.getArraySize());

    int fields = 0;
    for (final Iterator<Field> iter = info.getFieldIterator(); iter.hasNext(); ) {
      ++fields;
    }
    assertEquals(0, fields);

    int prims = 0;
    for (final Iterator<Field> iter = info.getNonprimitiveIterator(); iter.hasNext(); ) {
      ++prims;
    }
    assertEquals(0, prims);

  }

  public void testbooleanArray() {
    primArrayTst(new boolean[0], MemoryUsage.BOOLEAN_ARRAY);
  }

  public void testbyteArray() {
    primArrayTst(new byte[0], MemoryUsage.BYTE_ARRAY);
  }

  public void testcharArray() {
    primArrayTst(new char[0], MemoryUsage.CHAR_ARRAY);
  }

  public void testshortArray() {
    primArrayTst(new short[0], MemoryUsage.SHORT_ARRAY);
  }

  public void testintArray() {
    primArrayTst(new int[0], MemoryUsage.INT_ARRAY);
  }

  public void testfloatArray() {
    primArrayTst(new float[0], MemoryUsage.FLOAT_ARRAY);
  }

  public void testlongArray() {
    primArrayTst(new long[0], MemoryUsage.LONG_ARRAY);
  }

  public void testdoubleArray() {
    primArrayTst(new double[0], MemoryUsage.DOUBLE_ARRAY);
  }

  private static void primArrayTst(final Object obj, final int size) {
    final ClassMemory.Info info = ClassMemory.getMemoryInfo(obj.getClass());
    assertTrue(!info.isArray());
    assertTrue(info.isPrimitiveArray());

    assertEquals(info + "", MemoryUsage.ARRAY_SIZE, info.getSize());

    assertEquals(size, info.getArraySize());

    int fields = 0;
    for (final Iterator<Field> iter = info.getFieldIterator(); iter.hasNext(); ) {
      iter.next();
      ++fields;
    }
    assertEquals(0, fields);

    int prims = 0;
    for (final Iterator<Field> iter = info.getNonprimitiveIterator(); iter.hasNext(); ) {
      iter.next();
      ++prims;
    }
    assertEquals(0, prims);

  }

  public void testInteger() {
    final Integer x = 0;
    final ClassMemory.Info info = ClassMemory.getMemoryInfo(x.getClass());
    assertTrue(!info.isArray());
    assertTrue(!info.isPrimitiveArray());
    assertTrue(info.isObject());

    final int size = MemoryUsage.roundUp(MemoryUsage.OBJECT_SIZE + MemoryUsage.INT_SIZE);
    assertEquals(info + "", size, info.getSize());

    final int arr = 0;
    assertEquals(arr, info.getArraySize());

    int fields = 0;
    for (final Iterator<Field> iter = info.getFieldIterator(); iter.hasNext(); ) {
      iter.next();
      ++fields;
    }
    assertEquals(1, fields);

    int prims = 0;
    for (final Iterator<Field> iter = info.getNonprimitiveIterator(); iter.hasNext(); ) {
      iter.next();
      ++prims;
    }
    assertEquals(0, prims);

  }


  public void testIntegerArray() {
    final Integer[] x = new Integer[0];
    final ClassMemory.Info info = ClassMemory.getMemoryInfo(x.getClass());
    assertTrue(info.isArray());
    assertTrue(!info.isPrimitiveArray());

    final int size = MemoryUsage.ARRAY_SIZE;
    assertEquals(info + "", size, info.getSize());

    final int arr = MemoryUsage.REF_SIZE;
    assertEquals(arr, info.getArraySize());

    int fields = 0;
    for (final Iterator<Field> iter = info.getFieldIterator(); iter.hasNext(); ) {
      iter.next();
      ++fields;
    }
    assertEquals(0, fields);

    int prims = 0;
    for (final Iterator<Field> iter = info.getNonprimitiveIterator(); iter.hasNext(); ) {
      iter.next();
      ++prims;
    }
    assertEquals(0, prims);

  }
}


