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
package com.rtg.util.diagnostic;

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.PrintStream;

import junit.framework.TestCase;

/**
 */
public class TimerTest extends TestCase {

  /**
   * Test method for {@link com.rtg.util.diagnostic.Timer}.
   */
  public final void testForce() {
    final Timer ti = new Timer("B_timer");
    ti.start();
    ti.stop();
    ti.start();
    ti.stop();
    ti.setTime(1000000000);
    final String s = ti.toString();
    assertTrue(s.contains("Timer B_timer      1.00  count 2       0.50 bytes read 0"));
  }

  /**
   * Test method for {@link com.rtg.util.diagnostic.Timer}.
   * @throws IOException if an I/O error occurs.
   */
  public final void testLog1() throws IOException {
    final ByteArrayOutputStream ba = new ByteArrayOutputStream();
    final PrintStream pr = new PrintStream(ba);
    Diagnostic.setLogStream(pr);
    final Timer ti = new Timer("A_timer");
    checkUninitialized(ti);
    ti.reset();
    checkUninitialized(ti);

    ti.start();
    checkRunning(ti);
    ti.stop();
    checkStopped(ti);

    final String s1 = ti.toString();
    assertTrue(s1.contains("Timer A_timer "));
    assertTrue(s1.contains(" count 1 "));

    ti.start();
    checkRunning(ti);
    ti.stop();
    checkStopped(ti);

    final String s2 = ti.toString();
    assertTrue(s2.contains("Timer A_timer "));
    assertTrue(s2.contains(" count 2 "));
    ti.log();
    ti.log("foo");
    final String s3 = ba.toString();
    //System.err.println(s3);
    assertTrue(s3.contains(" Timer A_timer "));
    assertTrue(s3.contains(" Timer A_timer_foo "));
    assertTrue(s3.contains(" count 2 "));
    pr.close();
    ba.close();

    ti.reset();
    Diagnostic.setLogStream();
  }

  /**
   * Test method for {@link com.rtg.util.diagnostic.Timer}.
   * Supply names on reset.
   * @throws IOException if an I/O error occurs.
   */
  public final void testLog2() throws IOException {
    final ByteArrayOutputStream ba = new ByteArrayOutputStream();
    final PrintStream pr = new PrintStream(ba);
    Diagnostic.setLogStream(pr);
    final Timer ti = new Timer("A_timer");
    checkUninitialized(ti);
    ti.reset("B_timer");
    assertTrue(ti.toString().contains("Timer B_timer empty"));
    checkUninitialized(ti);

    ti.start();
    checkRunning(ti);
    ti.stop();
    checkStopped(ti);

    final String s1 = ti.toString();
    assertTrue(s1.contains("Timer B_timer "));
    assertTrue(s1.contains(" count 1 "));

    ti.start();
    checkRunning(ti);
    ti.stop(5);
    checkStopped(ti);
    ti.start();
    ti.stop(10);
    final String s2 = ti.toString();
    assertTrue(s2.contains("Timer B_timer "));
    assertTrue(s2.contains(" count 3 "));
    assertTrue(s2.contains("bytes read 15"));
    ti.log();
    ti.log("foo");
    final String s3 = ba.toString();
    //System.err.println(s3);
    assertTrue(s3.contains(" Timer B_timer "));
    assertTrue(s3.contains(" Timer B_timer_foo "));
    assertTrue(s3.contains(" count 3 "));
    pr.close();
    ba.close();

    ti.reset("C_timer");
    Diagnostic.setLogStream();
  }

  public final void testIncrement() {
    final Timer ti = new Timer("A_timer");
    checkUninitialized(ti);
    ti.increment(4200000000L);
    ti.increment(4200000000L);
    assertEquals("Timer A_timer      8.40  count 2       4.20 bytes read 0", ti.toString());
  }

  private void checkUninitialized(final Timer ti) {
    ti.integrity();
    try {
      ti.stop();
      fail();
    } catch (final IllegalStateException e) {
      //expected
    }
  }

  private void checkStopped(final Timer ti) {
    ti.integrity();
    try {
      ti.stop();
      fail();
    } catch (final IllegalStateException e) {
      //expected
    }
  }


  private void checkRunning(final Timer ti) {
    ti.integrity();
    try {
      ti.start();
      fail();
    } catch (final IllegalStateException e) {
      //expected
    }
    try {
      ti.toString();
      fail();
    } catch (final IllegalStateException e) {
      //expected
    }
  }

  public void testBadName() {
    try {
      new Timer("a b c");
      fail();
    } catch (final IllegalArgumentException e) {
      assertTrue(e.getMessage().startsWith("Name contains spaces:"));
    }
  }

  public void testEnum() {
    assertEquals("UNINITIALIZED", Timer.State.UNINITIALIZED.toString());
    //
    assertEquals(0, Timer.State.UNINITIALIZED.ordinal());
    assertEquals(1, Timer.State.STOPPED.ordinal());
    assertEquals(2, Timer.State.RUNNING.ordinal());
    assertEquals(Timer.State.UNINITIALIZED, Timer.State.valueOf("UNINITIALIZED"));
    assertEquals(Timer.State.STOPPED, Timer.State.valueOf("STOPPED"));
    assertEquals(Timer.State.RUNNING, Timer.State.valueOf("RUNNING"));

  }
}

