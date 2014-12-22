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
package com.rtg.util.cli;

import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;

/**
 * Tests the corresponding class.
 *
 */
public class WrappingStringBuilderTest extends TestCase {

  /**
   */
  public WrappingStringBuilderTest(final String name) {
    super(name);
  }

  public static Test suite() {
    return new TestSuite(WrappingStringBuilderTest.class);
  }

  /**
   * Main to run from tests from command line.
   * @param args ignored.
   */
  public static void main(final String[] args) {
    junit.textui.TestRunner.run(new TestSuite(WrappingStringBuilderTest.class));
  }

  private static final String LS = System.lineSeparator();

  public void test() {
    final WrappingStringBuilder b = new WrappingStringBuilder();
    assertEquals("", b.toString());
    assertEquals(b, b.append("hi"));
    assertEquals("hi", b.toString());
    assertEquals(b, b.append(" there peasant"));
    assertEquals("hi there peasant", b.toString());
    try {
      b.setWrapWidth(-1);
      fail();
    } catch (final IllegalArgumentException e) {
      assertEquals("Wrap width must be positive.", e.getMessage());
    }
    b.setWrapWidth(10);
    assertEquals("hi there peasant", b.toString());
    assertEquals("hi there peasant a a a a a a a a  a a aa a a a a a a a a a a a a a a a aa a a a a a a aa a a a a a a a a a a a a a a a a a aa  aa a a a a aaa a a a a aa a a a a a a a aa  ", b.wrapText(" a a a a a a a a  a a aa a a a a a a a a a a a a a a a aa a a a a a a aa a a a a a a a a a a a a a a a a a aa  aa a a a a aaa a a a a aa a a a a a a a aa  ").toString());
    assertEquals(b, b.append('z'));
    assertEquals("hi there peasant a a a a a a a a  a a aa a a a a a a a a a a a a a a a aa a a a a a a aa a a a a a a a a a a a a a a a a a aa  aa a a a a aaa a a a a aa a a a a a a a aa  z", b.toString());
  }

  public void testInitial() {
    assertEquals("hi", new WrappingStringBuilder("hi").toString());
    assertEquals(" a a a a a a a a  a a aa a a a a a a a a a a a a a a a aa a a a a a a aa a a a a a a a a a a a a a a a a a aa  aa a a a a aaa a a a a aa a a a a a a a aa  ", new WrappingStringBuilder(" a a a a a a a a  a a aa a a a a a a a a a a a a a a a aa a a a a a a aa a a a a a a a a a a a a a a a a a aa  aa a a a a aaa a a a a aa a a a a a a a aa  ").toString());
  }

  public void testWrapText() {
    final WrappingStringBuilder b = new WrappingStringBuilder();
    b.setWrapWidth(70);
    try {
      b.wrapText("hello\nthere");
      fail();
    } catch (final IllegalArgumentException e) {
      assertEquals("Input string cannot contain line breaks.", e.getMessage());
    }
    assertEquals(b, b.wrapText("hello"));
    final String s = b.wrapText(" a a a a a a a a  a a aa a a a a a a a a a a a a a a a aa a a a a a a aa a a a a a a a a a a a a a a a a a aa  aa a a a a aaa a a a a aa a a a a a a a aa  ").toString();
    assertEquals(s, "hello a a a a a a a a  a a aa a a a a a a a a a a a a a a a aa a a a" + LS + "a a a aa a a a a a a a a a a a a a a a a a aa  aa a a a a aaa a a a a" + LS + "aa a a a a a a a aa  ", s);
  }

  public void testWrapText2() {
    final WrappingStringBuilder b = new WrappingStringBuilder();
    b.setWrapWidth(70);
    b.setWrapIndent("POX");
    try {
      b.wrapText("hello\nthere");
      fail();
    } catch (final IllegalArgumentException e) {
      // ok
    }
    assertEquals(b, b.wrapText("hello"));
    assertEquals(b, b.wrapText(" a a a a a a a a  a a aa a a a a a a a a a a a a a a a aa a a a a a a aa a a a a a a a a a a a a a a a a a aa  aa a a a a aaa a a a a aa a a a a a a a aa  "));
    final String s = b.toString();
    assertEquals(s, "hello a a a a a a a a  a a aa a a a a a a a a a a a a a a a aa a a a" + LS + "POXa a a aa a a a a a a a a a a a a a a a a a aa  aa a a a a aaa a a" + LS + "POXa a aa a a a a a a a aa  ", s);
  }

  public void testWrapText3() {
    final WrappingStringBuilder b = new WrappingStringBuilder();
    b.setWrapWidth(70);
    b.setWrapIndent(5);
    try {
      b.wrapText("hello\nthere");
      fail();
    } catch (final IllegalArgumentException e) {
      // ok
    }
    assertEquals(b, b.wrapText("hello"));
    final String s = b.wrapText(" a a a a a a a a  a a aa a a a a a a a a a a a a a a a aa a a a a a a aa a a a a a a a a a a a a a a a a a aa  aa a a a a aaa a a a a aa a a a a a a a aa  ").toString();
    assertEquals(s, "hello a a a a a a a a  a a aa a a a a a a a a a a a a a a a aa a a a" + LS + "     a a a aa a a a a a a a a a a a a a a a a a aa  aa a a a a aaa a" + LS + "     a a a aa a a a a a a a aa  ", s);
  }

  public void testBorderlinePrefix() {
    final WrappingStringBuilder b = new WrappingStringBuilder();
    assertEquals(b, b.wrapText(""));
    assertEquals(b, b.wrapWord(""));
    b.setWrapWidth(40);
    b.setWrapIndent(20);
    assertEquals(b, b.wrapText("hello there my friend, squeamish ossifrage"));
    assertEquals("hello there my friend, squeamish" + LS + "                    ossifrage", b.toString());
    assertEquals(b, b.wrapText(""));
    assertEquals(b, b.wrapWord(""));
    assertEquals("hello there my friend, squeamish" + LS + "                    ossifrage", b.toString());
  }

  public void testBorderlinePrefix2() {
    final WrappingStringBuilder b = new WrappingStringBuilder();
    assertEquals(b, b.wrapText(""));
    assertEquals(b, b.wrapWord(""));
    b.setWrapWidth(40);
    b.setWrapIndent(21);
    assertEquals(b, b.wrapText("hello there my friend, squeamish ossifrage"));
    assertEquals(b, b.wrapText(""));
    assertEquals("hello there my friend, squeamish ossifrage", b.toString());
  }

  public void testBorderlinePrefix3() {
    final WrappingStringBuilder b = new WrappingStringBuilder();
    assertEquals(b, b.wrapText(""));
    assertEquals(b, b.wrapWord(""));
    b.setWrapWidth(40);
    b.setWrapIndent(21);
    assertEquals(b, b.wrapText("hello there my friend, squeamish"));
    assertEquals(b, b.wrapWord(" ossifrage"));
    assertEquals("hello there my friend, squeamish ossifrage", b.toString());
    assertEquals(b, b.wrapText(""));
    assertEquals("hello there my friend, squeamish ossifrage", b.toString());
  }

  public void testBorderlinePrefix4() {
    final WrappingStringBuilder b = new WrappingStringBuilder();
    assertEquals(b, b.wrapText(""));
    assertEquals(b, b.wrapWord(""));
    b.setWrapWidth(40);
    b.setWrapIndent(20);
    assertEquals(b, b.wrapText("hello there my friend, squeamish"));
    assertEquals(b, b.wrapWord(" ossifrage"));
    assertEquals("hello there my friend, squeamish" + LS + "                    ossifrage", b.toString());
  }

  public void testBorderlinePrefix5() {
    final WrappingStringBuilder b = new WrappingStringBuilder();
    assertEquals(b, b.wrapText(""));
    assertEquals(b, b.wrapWord(""));
    b.setWrapWidth(40);
    b.setWrapIndent(20);
    assertEquals(b, b.wrapText("hello there my friend, squeamish"));
    assertEquals(b, b.wrapWord("ossifrage"));
    assertEquals("hello there my friend, squeamish" + LS + "                    ossifrage", b.toString());
  }

  public void testException() {
    final WrappingStringBuilder b = new WrappingStringBuilder();
    try {
      b.wrapText("hello\nthere");
      fail();
    } catch (final IllegalArgumentException e) {
      // ok
    }
  }

  public void testEnd0() {
    final WrappingStringBuilder b = new WrappingStringBuilder();
    b.setWrapWidth(40);
    b.wrapText("0123456789012345678901234567890123456789");
    b.wrapText("a");
    assertEquals("0123456789012345678901234567890123456789" + LS + "a", b.toString());
  }

}
