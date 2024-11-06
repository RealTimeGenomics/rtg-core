/*
 * Copyright (c) 2018. Real Time Genomics Limited.
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the
 *    distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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
    final CreateParams c = new CreateParams(20, 33, 33, 31, false, true, false, false);
    assertEquals(322, IndexUtils.bytes(c));
  }

  private static final String EXPECTED_MEM_STR = ""
    + "\tMemory\tHash\t160" + StringUtils.LS
    + "\tMemory\tValue\t80" + StringUtils.LS
    + "\tMemory\tInitial_position\t18" + StringUtils.LS
    + "\tMemory\tBit_vector\t64" + StringUtils.LS
    + "";
  public void testMemoryTostring() {
    final CreateParams c = new CreateParams(20, 33, 33, 31, false, true, false, false);
    assertEquals(EXPECTED_MEM_STR, IndexUtils.memToString(c));
  }

  private static final String EXPECTED_MEMSTR = ""
  + "Memory Usage\tbytes\tlength" + StringUtils.LS
  + "\t\t160\t20\tHash" + StringUtils.LS
  + "\t\t80\t20\tValue" + StringUtils.LS
  + "\t\t18\t18\tInitial Position" + StringUtils.LS
  + "\t\t64\t512\tBit vector" + StringUtils.LS
  + "\t\t322\t\tTotal bytes" + StringUtils.LS;

  public void testMemoryString() {
    final CreateParams c = new CreateParams(20, 33, 33, 31, false, true, false, false);
    assertEquals(EXPECTED_MEMSTR, IndexUtils.memString(c));
  }
}
