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


import java.lang.ref.SoftReference;
import java.util.Stack;

import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;

/**
 * JUnit test class. <p>
 *
 */
public final class ObjectWalkerTest extends TestCase {

  static final class WalkableObject {
    WalkableObject() { } //prevent instantiation
    int mInt = 1;
    SoftReference<?> mRef = new SoftReference<>(new Object());
   // private WalkableObject mChild = null;
  }

  int mCounter = 0;

  WalkableObject mWalkable;

  /**
   * Constructor.
   *
   * @param s The test name.
   */
  public ObjectWalkerTest(final String s) {
    super(s);
  }

  public void testWalking() {
    final WalkableObject wo = new WalkableObject();
    final ObjectWalker ow = new ObjectWalkerImpl();

    // Simple walking.
    mCounter = 0;
    ow.walk(wo);
    int expected = 5;
    final String vmName = System.getProperty("java.vm.name", "");
    if (vmName.startsWith("IBM J9") || vmName.startsWith("BEA JRockit")) {
      expected = 3;
    }
    assertEquals(expected, mCounter);
  }

  public void testWalking2() {
    final WalkableObject wo = new WalkableObject();
    final ObjectWalker ow = new ObjectWalkerImpl();
    // Walking without reference following
    ow.setFollowReferences(false);
    mCounter = 0;
    ow.walk(wo);
    assertEquals(2, mCounter); // Object plus the reference
  }

  public void testWalking3() {
    final WalkableObject wo = new WalkableObject();
    final ObjectWalker ow = new ObjectWalkerImpl();
    ow.setFollowReferences(false);
    // Walking this test
    mCounter = 0;
    mWalkable = wo;
    assertEquals(1, mWalkable.mInt);
    assertNotNull(mWalkable.mRef);
   // assertNull(mWalkable.mChild);
    ow.walk(this);
    final int expected = 5;
    assertEquals(expected, mCounter); // Test bits plus the walkable.
  }

  public void testWalking4() {
    final WalkableObject wo = new WalkableObject();
    final ObjectWalker ow = new ObjectWalkerImpl();
    ow.setFollowReferences(false);
    // Walking with filter
    mCounter = 0;
    mWalkable = wo;
    ow.setWalkFilter(new com.rtg.util.memory.ExcludeClassFilter(new Class<?>[0]) {
        @Override
        public boolean accept(final Stack<Object> path, final Object value) {
          return path.isEmpty() || !checkIsTest(path.peek());
        }
      });
    ow.walk(this);
    final int expected = 3;
    assertEquals(expected, mCounter); // Stop walking past the test
  }

  static boolean checkIsTest(final Object o) {
    return o instanceof Test;
  }

  @Override
  public void tearDown() {
    mWalkable = null;
  }

  /**
   * Adds tests to suite to be run by main
   *
   * @return The test suite.
   */
  public static Test suite() {
    return new TestSuite(ObjectWalkerTest.class);
  }


  /**
   * Main method needed to make a self runnable class
   *
   * @param args The command line arguments.
   */
  /**
   * Main to run from tests from command line.
   * @param args ignored.
   */
  public static void main(final String[] args) {
    junit.textui.TestRunner.run(suite());
  }

  private class ObjectWalkerImpl extends ObjectWalker {

    public ObjectWalkerImpl() {
    }

    @Override
    protected void visitObject(final Stack<Object> path, final Object value) {
      // super.visitObject(path, value);  // Uncomment for trace on system err
      mCounter++;
    }
  }
}


