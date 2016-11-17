/*
 * Copyright (c) 2016. Real Time Genomics Limited.
 *
 * Use of this source code is bound by the Real Time Genomics Limited Software Licence Agreement
 * for Academic Non-commercial Research Purposes only.
 *
 * If you did not receive a license accompanying this file, a copy must first be obtained by email
 * from support@realtimegenomics.com.  On downloading, using and/or continuing to use this source
 * code you accept the terms of that license agreement and any amendments to those terms that may
 * be made from time to time by Real Time Genomics Limited.
 */

package com.rtg.variant.avr;

import com.rtg.util.TestUtils;
import com.rtg.vcf.header.MetaType;

import junit.framework.TestCase;

/**
 */
public class AnnotationDataTypeTest extends TestCase {

  public void testEnum() {
    TestUtils.testEnum(AnnotationDataType.class, "[BOOLEAN, INTEGER, DOUBLE, STRING]");
    assertEquals(0, AnnotationDataType.BOOLEAN.ordinal());
    assertEquals(1, AnnotationDataType.INTEGER.ordinal());
    assertEquals(2, AnnotationDataType.DOUBLE.ordinal());
    assertEquals(3, AnnotationDataType.STRING.ordinal());

    assertEquals(4, AnnotationDataType.values().length);
  }

  public void testGetClassType() {
    assertEquals(String.class, AnnotationDataType.STRING.getClassType());
    assertEquals(Double.class, AnnotationDataType.DOUBLE.getClassType());
    assertEquals(Boolean.class, AnnotationDataType.BOOLEAN.getClassType());
    assertEquals(Integer.class, AnnotationDataType.INTEGER.getClassType());
  }

  public void testIsMetaTypeCompatiable() {
    assertTrue(AnnotationDataType.STRING.isMetaTypeCompatible(MetaType.STRING));
    assertTrue(AnnotationDataType.STRING.isMetaTypeCompatible(MetaType.CHARACTER));
    assertFalse(AnnotationDataType.STRING.isMetaTypeCompatible(MetaType.INTEGER));
    assertFalse(AnnotationDataType.STRING.isMetaTypeCompatible(MetaType.FLAG));
    assertFalse(AnnotationDataType.STRING.isMetaTypeCompatible(MetaType.FLOAT));

    assertTrue(AnnotationDataType.INTEGER.isMetaTypeCompatible(MetaType.INTEGER));
    assertFalse(AnnotationDataType.INTEGER.isMetaTypeCompatible(MetaType.CHARACTER));
    assertFalse(AnnotationDataType.INTEGER.isMetaTypeCompatible(MetaType.STRING));
    assertFalse(AnnotationDataType.INTEGER.isMetaTypeCompatible(MetaType.FLAG));
    assertFalse(AnnotationDataType.INTEGER.isMetaTypeCompatible(MetaType.FLOAT));

    assertTrue(AnnotationDataType.DOUBLE.isMetaTypeCompatible(MetaType.FLOAT));
    assertFalse(AnnotationDataType.DOUBLE.isMetaTypeCompatible(MetaType.CHARACTER));
    assertFalse(AnnotationDataType.DOUBLE.isMetaTypeCompatible(MetaType.INTEGER));
    assertFalse(AnnotationDataType.DOUBLE.isMetaTypeCompatible(MetaType.FLAG));
    assertFalse(AnnotationDataType.DOUBLE.isMetaTypeCompatible(MetaType.STRING));

    assertTrue(AnnotationDataType.BOOLEAN.isMetaTypeCompatible(MetaType.FLAG));
    assertFalse(AnnotationDataType.BOOLEAN.isMetaTypeCompatible(MetaType.CHARACTER));
    assertFalse(AnnotationDataType.BOOLEAN.isMetaTypeCompatible(MetaType.INTEGER));
    assertFalse(AnnotationDataType.BOOLEAN.isMetaTypeCompatible(MetaType.STRING));
    assertFalse(AnnotationDataType.BOOLEAN.isMetaTypeCompatible(MetaType.FLOAT));
  }

  public void testConvert() {
    for (AnnotationDataType adt : AnnotationDataType.values()) {
      assertNull(adt.stringToObjectOfType(null));
      if (adt != AnnotationDataType.STRING && adt != AnnotationDataType.BOOLEAN) {
        try {
          adt.stringToObjectOfType("adsafsfd");
          fail(adt + " formatted bad value");
        } catch (IllegalArgumentException iae) {
          // expected
        }
      }
    }
    Object o = AnnotationDataType.STRING.stringToObjectOfType("banana");
    assertTrue(o instanceof String);
    assertEquals("banana", o);

    o = AnnotationDataType.INTEGER.stringToObjectOfType("123");
    assertTrue(o instanceof Integer);
    assertEquals(123, o);

    o = AnnotationDataType.DOUBLE.stringToObjectOfType("123");
    assertTrue(o instanceof Double);
    assertEquals(123.0, o);

    o = AnnotationDataType.BOOLEAN.stringToObjectOfType("true");
    assertTrue(o instanceof Boolean);
    assertEquals(Boolean.TRUE, o);

  }
}
