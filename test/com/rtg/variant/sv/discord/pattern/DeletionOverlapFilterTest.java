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

package com.rtg.variant.sv.discord.pattern;

import com.rtg.variant.sv.discord.DiscordBedRecord;

import junit.framework.TestCase;

/**
 *         Date: 20/03/12
 *         Time: 1:26 PM
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
    DeletionOverlapFilter filter = new DeletionOverlapFilter();
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
      i++;
    }
    assertEquals(EXPECTED.length, i);
  }
}
