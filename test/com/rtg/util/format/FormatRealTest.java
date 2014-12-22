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
package com.rtg.util.format;


import junit.framework.TestCase;

/**
 * Tests FormatReal. Run from the command line with:<p>
 *
 * <code>
 * java junit.textui.TestRunner com.rtg.util.format.FormatRealTest<br>
 * java junit.swingui.TestRunner com.rtg.util.format.FormatRealTest<br>
 * java com.rtg.util.format.FormatRealTest<br>
 * </code>
 *
 */
public class FormatRealTest extends TestCase {

  public FormatRealTest(final String name) {
    super(name);
  }

  @Override
  public void setUp() {
    clearStrBuffer();
  }

  @Override
  public void tearDown() {
    mStrBuff = null;
    mFReal = null;
  }


  /** Clear StringBuilder object for test methods  */
  private void clearStrBuffer() {
    mStrBuff = new StringBuilder();
  }


  /**
   * Tests constructor and format methods FormatReal constructor is
   * sent # of places to left of decimal and # of spaces to right of
   * decimal as integer parameters
   */
  public void testFormatDoubles() {
    mFReal = new FormatReal(2, 3);
    mStr = mFReal.format(mStrBuff, 14.7).toString();
    assertEquals("Should be equal", "14.700", mStr);

    mFReal = new FormatReal(2, 1);
    clearStrBuffer();
    mStr = mFReal.format(mStrBuff, -4.75).toString();
    // Should round decimal up to 1st place
    assertEquals("Should be equal", "-4.8", mStr);

    mFReal = new FormatReal(1, 2);
    clearStrBuffer();
    mStr = mFReal.format(mStrBuff, -23.27).toString();
    assertEquals("Should be equal", "#-23.27#", mStr);

    mFReal = new FormatReal(2, 2);
    clearStrBuffer();
    mStr = mFReal.format(mStrBuff, Float.MIN_VALUE).toString();
    assertEquals("Should be equal", "0.00", mStr.trim()); // it will round MIN to zero

    mFReal = new FormatReal(2, 5);
    clearStrBuffer();
    mStr = mFReal.format(mStrBuff, Float.MAX_VALUE).toString();
    //System.err.println("MAX >>> " + mStr + mStr.substring(mStr.indexOf(".")+1 ).length());
    assertEquals("Should be the size of 6 ", 6, mStr.substring(mStr.indexOf(".") + 1).length());


  }


  /**
   * Tests constructor and format methods FormatReal constructor is
   * sent # of places to left of decimal and # of spaces to right of
   * decimal as integer parameters
   */
  public void testFormatFloats() {
    mFReal = new FormatReal(2, 3);
    mStr = mFReal.format(mStrBuff, 14.7f).toString();
    assertEquals("Should be equal", "14.700", mStr);

    mFReal = new FormatReal(1, 1);
    clearStrBuffer();
    mStr = mFReal.format(mStrBuff, 4.75f).toString();
    // Should round decimal up to 1st place
    assertEquals("Should be equal", "4.8", mStr);

    mFReal = new FormatReal(2, 10);
    clearStrBuffer();
    mStr = mFReal.format(mStrBuff, -3.2e10f).toString();
    assertEquals("Should be equal", "#-32000000000.0000000000#", mStr);

    mFReal = new FormatReal(2, 2);
    clearStrBuffer();
    mStr = mFReal.format(mStrBuff, Float.MIN_VALUE).toString();
    assertEquals("Should be equal", " 0.00", mStr); // it will round MIN to zero


    mFReal = new FormatReal(5, 5);
    clearStrBuffer();
    mStr = mFReal.format(mStrBuff, Float.MAX_VALUE).toString();

    //System.err.println("MAX <<< " + mStr + str.substring(mStr.indexOf(".")+1 ).length());
    assertEquals("Should be the size of 6 ", 6, mStr.substring(mStr.indexOf(".") + 1).length());

  }


  /** Sends negative values to constructor  */
  public void testInvalidConstructor() {
    // Negative # of places to left of decimal specified
    try {
      mFReal = new FormatReal(-1, 1);
      fail("should throw IllegalArgumentException");
    } catch (final IllegalArgumentException e) { }
  }


  /** Test with 0 decimal places  */
  public void testConstructorZeroDecimalPlaces() {
    mFReal = new FormatReal(2, 0);
    mStr = mFReal.format(mStrBuff, 14.7).toString();
    assertEquals("Should be equal", "15.", mStr);
  }


  /**
   * Test with the negative input expected : uncaught
   * NegativeArraySizeException
   */
  public void testConstructorNegativeInput() {
    try {
      mFReal = new FormatReal(-100, 0);
      fail("should throw illegal argument exception");
    } catch (final IllegalArgumentException e) { }
  }



  /**
   * Tests format method with null input expected uncaught
   * NullPointerException
   */
  public void testFormatNull() {
    mFReal = new FormatReal(2, 2);
    try {
      mStrBuff = mFReal.format(null, 0.0d);
      fail("should throw exception");
    } catch (final RuntimeException e) { }
  }


  /** Tests format method with NaN  */
  public void testFormatNaN() {
    mFReal = new FormatReal(2, 2);
    mStr = mFReal.format(mStrBuff, Float.NaN).toString();
    assertEquals("Should be NaN", mStr, "  NaN");
    mStrBuff = new StringBuilder();
    mStr = mFReal.format(mStrBuff, Float.NaN).toString();
    assertEquals("Should be NaN", mStr, "  NaN");
  }



  /**
   * Tests format method with Pos/Neg Infinities Is Infinity argument
   * legal?
   */
  public void testFormatInfinityDP() {
    mFReal = new FormatReal(2, 2);
    // Test using doubles
    mStr = mFReal.format(mStrBuff, Float.POSITIVE_INFINITY).toString();
    assertEquals("Should be Infinity", "#Infinity#", mStr);
  }


  /**
   * Tests format method with Pos/Neg Infinities Is Infinity argument
   * legal?
   */
  public void testFormatInfinityDN() {
    mFReal = new FormatReal(2, 2);
    clearStrBuffer();
    mStr = mFReal.format(mStrBuff, Float.NEGATIVE_INFINITY).toString();
    assertEquals("Should be Infinity", "#-Infinity#", mStr);
  }


  /**
   * Tests format method with Pos/Neg Infinities Is Infinity argument
   * legal?
   */
  public void testFormatInfinitySP() {
    mFReal = new FormatReal(2, 2);
    clearStrBuffer();
    mStr = mFReal.format(mStrBuff, Float.POSITIVE_INFINITY).toString();

    assertEquals("Should be Infinity", "#Infinity#", mStr);
  }



  /**
   * Tests format method with Pos/Neg Infinities Is Infinity argument
   * legal?
   */
  public void testFormatInfinitySN() {
    mFReal = new FormatReal(2, 2);
    clearStrBuffer();
    mStr = mFReal.format(mStrBuff, Float.NEGATIVE_INFINITY).toString();

    assertEquals("Should be negative Infinity", "#-Infinity#", mStr);
  }

  private StringBuilder mStrBuff;
  private String mStr;
  private FormatReal mFReal;
}
