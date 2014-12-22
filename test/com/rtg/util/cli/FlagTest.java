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

import java.util.Arrays;

import com.rtg.util.EnumHelper;
import com.rtg.util.PseudoEnum;

import junit.framework.TestCase;

/**
 * Test Flag
 */
public class FlagTest extends TestCase {

  static class DummyEnum implements PseudoEnum {
    public static final DummyEnum WORD = new DummyEnum(0, "WORD");
    public static final DummyEnum SEGMENT = new DummyEnum(1, "SEGMENT");

    private final int mOrdinal;
    private final String mName;

    public DummyEnum(final int ordinal, final String name) {
      mOrdinal = ordinal;
      mName = name;
    }
    @Override
    public String name() {
      return mName;
    }
    @Override
    public int ordinal() {
      return mOrdinal;
    }
    @Override
    public String toString() {
      return mName;
    }
    private static final EnumHelper<DummyEnum> HELPER = new EnumHelper<>(DummyEnum.class, new DummyEnum[] {WORD, SEGMENT});

    public static DummyEnum[] values() {
      return HELPER.values();
    }

    public static DummyEnum valueOf(final String str) {
      return HELPER.valueOf(str);
    }
  }

  public void testIsValidEnum() {
    assertTrue(Flag.isValidEnum(DummyEnum.class));

    assertFalse(Flag.isValidEnum(String.class));
    assertFalse(Flag.isValidEnum(null));
  }

  public void testValues2() {
    final String[] values = Flag.values(DummyEnum.class);
    assertNotNull(values);
    assertEquals("[word, segment]", Arrays.toString(values));
    assertEquals(DummyEnum.WORD, Flag.valueOf(DummyEnum.class, "WORD"));
    assertEquals(DummyEnum.WORD, Flag.instanceHelper(DummyEnum.class, "WORD"));
    assertEquals(DummyEnum.WORD, Flag.instanceHelper(DummyEnum.class, "word"));
    assertEquals(DummyEnum.WORD, Flag.instanceHelper(DummyEnum.class, "WorD"));
  }

  public void testValues3() {
    assertNull(Flag.values(String.class));
    assertNull(Flag.values(null));
  }

  public void testHashCode() {
    Flag anon = new Flag(null, null, "anonymous flag", 1, 1, Integer.class, "int", null, "");
    assertEquals(0, anon.hashCode());
  }

  public void testMinMax() {
    assertEquals("", Flag.minMaxUsage(0, 1)); //this will be dealt with by the normal optional usage
    assertEquals("", Flag.minMaxUsage(1, 1)); //this will be dealt with by the normal required usage
    assertEquals("May be specified up to 2 times", Flag.minMaxUsage(0, 2));
    assertEquals("May be specified up to 3 times", Flag.minMaxUsage(0, 3));
    assertEquals("May be specified 0 or more times", Flag.minMaxUsage(0, Integer.MAX_VALUE));
    assertEquals("Must be specified 1 or 2 times", Flag.minMaxUsage(1, 2));
    assertEquals("Must be specified 1 to 3 times", Flag.minMaxUsage(1, 3));
    assertEquals("Must be specified 1 or more times", Flag.minMaxUsage(1, Integer.MAX_VALUE));
    assertEquals("Must be specified 2 times", Flag.minMaxUsage(2, 2));
    assertEquals("Must be specified 2 or 3 times", Flag.minMaxUsage(2, 3));
    assertEquals("Must be specified 2 to 4 times", Flag.minMaxUsage(2, 4));
    assertEquals("Must be specified 2 or more times", Flag.minMaxUsage(2, Integer.MAX_VALUE));
  }
}
