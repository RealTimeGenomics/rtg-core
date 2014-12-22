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
package com.rtg.ml;

import junit.framework.TestCase;

/**
 */
public class AttributeTest extends TestCase {

  public AttributeTest(final String name) {
    super(name);
  }

  public void test() {
    final Attribute a = new Attribute("higgs", MlDataType.DOUBLE);
    assertEquals("higgs", a.getName());
    assertEquals(MlDataType.DOUBLE, a.getDataType());
  }

  public void testEncodingDouble() {
    final Attribute d = new Attribute("tisadouble", MlDataType.DOUBLE);
    final double orig = 3.0;
    double vala = d.encodeValue(orig);
    double valb = d.encodeValue(orig);
    assertEquals(vala, valb);
    assertEquals(orig, d.decodeValue(vala));
    assertEquals(orig, d.decodeValue(valb));
  }

  public void testEncodingBoolean() {
    final Attribute d = new Attribute("tisabool", MlDataType.BOOLEAN);
    final boolean orig = true;
    double vala = d.encodeValue(orig);
    double valb = d.encodeValue(orig);
    assertEquals(vala, valb);
    assertEquals(orig, d.decodeValue(vala));
    assertEquals(orig, d.decodeValue(valb));
  }

  public void testEncodingInteger() {
    final Attribute d = new Attribute("tisanint", MlDataType.INTEGER);
    final int orig = 555356;
    double vala = d.encodeValue(orig);
    double valb = d.encodeValue(orig);
    assertEquals(vala, valb);
    assertEquals(orig, d.decodeValue(vala));
    assertEquals(orig, d.decodeValue(valb));
  }

  public void testEncodingString() {
    final Attribute d = new Attribute("tisastring", MlDataType.STRING);
    //load up a bunch of values
    d.encodeValue("stringa");
    d.encodeValue("stringb");
    d.encodeValue("stringc");
    d.encodeValue("stringd");
    d.encodeValue("stringe");
    final String orig = "stringc";
    double vala = d.encodeValue(orig);
    double valb = d.encodeValue(orig);
    assertEquals(vala, valb);
    assertEquals(orig, d.decodeValue(vala));
    assertEquals(orig, d.decodeValue(valb));
  }

}
