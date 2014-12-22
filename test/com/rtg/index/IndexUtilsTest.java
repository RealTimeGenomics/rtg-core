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
package com.rtg.index;

import com.rtg.index.params.CreateParams;
import com.rtg.util.StringUtils;

import junit.framework.TestCase;

/**
 *
 *
 *
 */
public class IndexUtilsTest extends TestCase {

  /**
   */
  public IndexUtilsTest(final String name) {
    super(name);
  }

  public void testBytes() {
    final CreateParams c = new CreateParams(20, 33, 33, false, false, false);
    assertEquals(322, IndexUtils.bytes(c));
  }

  private static final String EXPECTED_MEM_STR = ""
    + "\tMemory\tHash\t160" + StringUtils.LS
    + "\tMemory\tValue\t80" + StringUtils.LS
    + "\tMemory\tInitial_position\t18" + StringUtils.LS
    + "\tMemory\tBit_vector\t64" + StringUtils.LS
    + "";
  public void testMemoryTostring() {
    final CreateParams c = new CreateParams(20, 33, 33, false, false, false);
    assertEquals(EXPECTED_MEM_STR, IndexUtils.memToString(c));
  }

  private static final String EXPECTED_MEMSTR = StringUtils.LS
  + "Memory Usage\tbytes\tlength" + StringUtils.LS
  + "\t\t160\t20\tHash" + StringUtils.LS
  + "\t\t80\t20\tValue" + StringUtils.LS
  + "\t\t18\t18\tInitial Position" + StringUtils.LS
  + "\t\t64\t512\tBit vector" + StringUtils.LS
  + "\t\t322\t\tTotal" + StringUtils.LS;

  public void testMemoryString() {
    final CreateParams c = new CreateParams(20, 33, 33, false, false, false);
    assertEquals(EXPECTED_MEMSTR, IndexUtils.memString(c));
  }
}
