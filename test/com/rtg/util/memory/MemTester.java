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


/**
 * Measures the memory used by primitives and some simple classes. Does
 * this by carefully allocating (and not deallocating) a sequence of
 * array. Must be run with memory upper and lower bounds set so garbage
 * collection never kicks in. Objects up to 2K seem to be added and
 * subtracted during this. java -Xms100m -Xmx100m
 * com.rtg.util.memory.MemTester [-v] The -v option prints a memory report
 * for each object that is tested.
 *
 */
//WARNING: I apologize for the use of global variables and other hacky stuff in here.
//Some of it is necessary to.create this work reliably.
//JGC.
public class MemTester {
  /** Number of repetitions of objects when measuring size */
  static final int LEN = 100000;


  static long sLastFree;


  static long sFree;


  protected static final Runtime RUNTIME = Runtime.getRuntime();

  /**
   * Fraction by which an estimate has to differ from actual measured
   * before an error is recorded - needed because not all memory usage
   * can be accounted for.
   */
  static final double TOLERANCE = 0.2; //Heuristic - may need to vary with JVM

  /**
   * Holds output which is written all at once at end to prevent
   * creating and destroying memory.
   */
  static final StringBuilder OUT = new StringBuilder(10000); //big enough to hold all output

  /**
   * Keeps track of which pass through the memory test - different
   * things are done on each pass: 0 dummy to let memory system settle
   * down 1 check free memory against estimated size 2 measure size
   * against exact memory size calculation
   */
  static int sPass = 0;

  /**
   * Amount of printing - 0 nothing unless errors - 1 each test - 2
   * each test plus memory information.
   */
  static int sVerbose = 0;

  /** Keeps count of the number of errors. */
  static int sErrors = 0;


  /**
   * Perform the check.
   *
   * @param msg message
   * @param check check value
   * @param obj object to check
   */
  public static void check(final String msg, final int check, final Object obj) {
    //first pass seems to give inconsistent numbers so do not test on it
    if (sPass == 0) {
      return;
    }
    final long d;
    if (sPass == 1) {
      //check against free mem
      sFree = RUNTIME.freeMemory();
      d = sLastFree - sFree - MemoryUsage.roundUp(MemoryUsage.ARRAY_SIZE);
    } else {
      //check against memSize
      final MemoryUsage mem = new MemoryUsage(obj);
      if (sVerbose == 2) {
        System.err.println(mem);
      }
      d = mem.getSize();
      final long s = MemoryUsage.size(obj);
      if (d != s) {
        throw new RuntimeException("simple version of memory size disagrees:"
          + d + ":" + s + "\n" + mem);
      }
    }

    final double estimate = d / (double) LEN;
    OUT.append(msg).append(estimate);
    if (check > estimate + TOLERANCE || check < estimate - TOLERANCE) {
      sErrors++;
      OUT.append(" disagrees with ").append(check).append(" Error **");
    }
    OUT.append("\n");
    sLastFree = sFree;
  }


  /**
   * Main progream.
   *
   * @param args verbosity
   */
  public static void main(final String[] args) {
    if (args.length != 0) {
      if (args[0].equals("-v")) {
        sVerbose = 1;
      } else if (args[0].equals("-vv")) {
        sVerbose = 2;
      }
    }
    OUT.append("starting memory measurements\n");
    final MemTester mem = new MemTester();
    sPass = 0;
    mem.memCheck();

    OUT.append("comparison with free memory calculation\n");
    sPass = 1; //check freemem
    mem.memCheck();
    OUT.append("\n");

    OUT.append("comparison with memSize calculation\n");
    sPass = 2; //check memsize
    mem.memCheck();

    if (sVerbose >= 1 || sErrors > 0) {
      System.out.print(OUT.toString());
    }
  }


  static final class Thing {

  }


  static class Thing1 {
    protected int mDummy;
    protected int mDummy2;
  }


  static class Thing2 {
    protected byte mDummy;
    protected byte mDummy2;
  }


  static class Nest0 {
    protected byte mDummy;
  }


  static class Nest1 extends Nest0 {
    protected byte mDummy1;
  }


  static class Nest2 extends Nest1 {
    protected byte mDummy2;
  }


  static class Nest3 extends Nest2 {
    protected byte mDummy3;
  }


  static class Nest4 extends Nest3 {
    protected byte mDummy4;
  }

  //Check how many bytes for a byte
  protected static class By {
    protected byte mD1, mD2, mD3, mD4, mD5, mD6, mD7, mD8;
  }

  //Check how many bytes for an int
  protected static class In {
    protected int mD1, mD2, mD3, mD4, mD5, mD6, mD7, mD8;
  }


  class Inner0 {
  }


  Inner1 mInner1 = new Inner1();


  static class Inner1 {
    protected int mI;


    class InnerInner {
      protected int mJ;
    }


    InnerInner[] create() {
      final InnerInner[] r = new InnerInner[LEN];
      for (int i = 0; i < LEN; i++) {
        r[i] = new InnerInner();
      }
      return r;
    }
  }


  class Inner2 {
    protected int mI0, mI1;
  }

  static void updateLastFree() {
    RUNTIME.gc();
    sLastFree = RUNTIME.freeMemory();
  }

  void memCheck() {
    updateLastFree();

    check("nothing:", 0, new Object[0]);
    check("boolean:", MemoryUsage.BOOLEAN_ARRAY, new boolean[LEN]);
    check("byte:", MemoryUsage.BYTE_ARRAY, new byte[LEN]);
    check("char:", MemoryUsage.CHAR_ARRAY, new char[LEN]);
    check("short:", MemoryUsage.SHORT_ARRAY, new short[LEN]);
    check("int:", MemoryUsage.INT_ARRAY, new int[LEN]);
    check("float:", MemoryUsage.FLOAT_ARRAY, new float[LEN]);
    check("long:", MemoryUsage.LONG_ARRAY, new long[LEN]);
    check("double:", MemoryUsage.DOUBLE_ARRAY, new double[LEN]);
    check("object reference:", MemoryUsage.REF_SIZE, new Object[LEN]);

    memCheckForLoops1();
    memCheckForLoops2();
  }

  private void memCheckForLoops1() {
    final Thing[] oo = new Thing[LEN];
    for (int i = 0; i < LEN; i++) {
      oo[i] = new Thing();
    }
    check("empty object:", roundUp(0), oo);

    final Thing1[] o1 = new Thing1[LEN];
    for (int i = 0; i < LEN; i++) {
      o1[i] = new Thing1();
    }
    check("object with two ints:", roundUp(2 * MemoryUsage.INT_SIZE), o1);

    final Thing2[] o2 = new Thing2[LEN];
    for (int i = 0; i < LEN; i++) {
      o2[i] = new Thing2();
    }
    check("object with two bytes:", roundUp(2 * MemoryUsage.BYTE_SIZE), o2);

    final Nest0[] n0 = new Nest0[LEN];
    for (int i = 0; i < LEN; i++) {
      n0[i] = new Nest0();
    }
    check("object Nest0:", roundUp(MemoryUsage.BYTE_SIZE), n0);

    final Nest1[] n1 = new Nest1[LEN];
    for (int i = 0; i < LEN; i++) {
      n1[i] = new Nest1();
    }
    check("object Nest1:", roundUp(2 * MemoryUsage.BYTE_SIZE), n1);

    final Nest2[] n2 = new Nest2[LEN];
    for (int i = 0; i < LEN; i++) {
      n2[i] = new Nest2();
    }
    check("object Nest2:", roundUp(3 * MemoryUsage.BYTE_SIZE), n2);

    final Nest3[] n3 = new Nest3[LEN];
    for (int i = 0; i < LEN; i++) {
      n3[i] = new Nest3();
    }
    check("object Nest3:", roundUp(4 * MemoryUsage.BYTE_SIZE), n3);

    final Nest4[] n4 = new Nest4[LEN];
    for (int i = 0; i < LEN; i++) {
      n4[i] = new Nest4();
    }
    check("object Nest4:", roundUp(5 * MemoryUsage.BYTE_SIZE), n4);

    final By[] bz = new By[LEN];
    for (int i = 0; i < LEN; i++) {
      bz[i] = new By();
    }
    check("8 bytes:", roundUp(8 * MemoryUsage.BYTE_SIZE), bz);

    final In[] in = new In[LEN];
    for (int i = 0; i < LEN; i++) {
      in[i] = new In();
    }
    check("8 ints:", roundUp(8 * MemoryUsage.INT_SIZE), in);
  }

  private void memCheckForLoops2() {
    final Object[][] a0 = new Object[LEN][];
    for (int i = 0; i < LEN; i++) {
      a0[i] = new Object[0];
    }
    check("object array:", ra(0), a0);

    final byte[][] b0 = new byte[LEN][];
    for (int i = 0; i < LEN; i++) {
      b0[i] = new byte[0];
    }
    check("byte array:", ra(0), b0);

    final byte[][] b1 = new byte[LEN][];
    for (int i = 0; i < LEN; i++) {
      b1[i] = new byte[1];
    }
    check("byte array[1]:", ra(MemoryUsage.BYTE_ARRAY), b1);

    final byte[][] b4 = new byte[LEN][];
    for (int i = 0; i < LEN; i++) {
      b4[i] = new byte[4];
    }
    check("byte array[4]:", ra(4 * MemoryUsage.BYTE_ARRAY), b4);

    final byte[][] b5 = new byte[LEN][];
    for (int i = 0; i < LEN; i++) {
      b5[i] = new byte[5];
    }
    check("byte array[5]:", ra(5 * MemoryUsage.BYTE_ARRAY), b5);

    final Inner0[] i0 = new Inner0[LEN];
    for (int i = 0; i < LEN; i++) {
      i0[i] = new Inner0();
    }
    check("inner0:", roundUp(MemoryUsage.INNER_SIZE), i0);

    final Inner1[] i1 = new Inner1[LEN];
    for (int i = 0; i < LEN; i++) {
      i1[i] = new Inner1();
    }
    check("inner1:", roundUp(MemoryUsage.INT_SIZE + MemoryUsage.INNER_SIZE), i1);

    final Inner2[] i2 = new Inner2[LEN];
    for (int i = 0; i < LEN; i++) {
      i2[i] = new Inner2();
    }
    check("inner2:", roundUp(2 * MemoryUsage.INT_SIZE + MemoryUsage.INNER_SIZE), i2);

    final Inner1.InnerInner[] ia = mInner1.create();
    check("inner inner:", roundUp(MemoryUsage.INT_SIZE + MemoryUsage.INNER_SIZE), ia);

    final Integer[] oi = new Integer[LEN];
    for (int i = 0; i < LEN; i++) {
      oi[i] = 0;
    }
    check("Integer:", roundUp(MemoryUsage.INT_SIZE), oi);
  }


  /**
   * Total size of an object including extra fields.
   *
   * @param size of extra fields.
   * @return total size.
   */
  private int roundUp(final int size) {
    return MemoryUsage.REF_SIZE + MemoryUsage.roundUp(MemoryUsage.OBJECT_SIZE + size);
  }


  /**
   * Total size of an array including dope vector.
   *
   * @param size of dope vector.
   * @return total size.
   */
  private int ra(final int size) {
    return MemoryUsage.REF_SIZE + MemoryUsage.roundUp(MemoryUsage.ARRAY_SIZE + size);
  }

}
