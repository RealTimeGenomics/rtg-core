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

package com.rtg.variant.sv.discord.pattern;

import com.rtg.variant.sv.discord.DiscordBedRecord;

import junit.framework.TestCase;

/**
 */
public class DeletionOverlapFilterTest extends TestCase {

  static final String[] EXPECTED = {
      "bar\t10\t30"
      , "foo\t12\t19"
      , "foo\t21\t29"
      , "fuz\t9\t13"
      , "goo\t5\t10"
      , "goo\t10\t20"
      , "heh\t20\t22"
      , "heh\t10\t20"
  };
  public void testOverlap() {
    final DeletionOverlapFilter filter = new DeletionOverlapFilter();
    filter.add(new DiscordBedRecord("foo", 10, 30));
    filter.add(new DiscordBedRecord("foo", 12, 19));
    filter.add(new DiscordBedRecord("bar", 10, 30));
    filter.add(new DiscordBedRecord("foo", 21, 29));
    filter.add(new DiscordBedRecord("foo", 12, 30));
    filter.add(new DiscordBedRecord("bar", 1, 25));

    filter.add(new DiscordBedRecord("fuz", 5, 10));
    filter.add(new DiscordBedRecord("fuz", 9, 13));
    filter.add(new DiscordBedRecord("fuz", 10, 20));

    filter.add(new DiscordBedRecord("goo", 5, 10));
    filter.add(new DiscordBedRecord("goo", 9, 20));
    filter.add(new DiscordBedRecord("goo", 10, 20));

    filter.add(new DiscordBedRecord("heh", 20, 22));
    filter.add(new DiscordBedRecord("heh", 10, 21));
    filter.add(new DiscordBedRecord("heh", 10, 20));
    int i = 0;
    for (DiscordBedRecord b : filter) {
      assertEquals(EXPECTED[i], b.toString());
//      System.err.println(EXPECTED[i] + " == " + b.toString());
      ++i;
    }
    assertEquals(EXPECTED.length, i);
  }
}
