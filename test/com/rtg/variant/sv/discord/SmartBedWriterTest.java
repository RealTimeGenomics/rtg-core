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

package com.rtg.variant.sv.discord;

import static com.rtg.util.StringUtils.LS;

import java.io.ByteArrayOutputStream;
import java.io.IOException;

import com.rtg.util.TestUtils;

import junit.framework.TestCase;


/**
 */
public class SmartBedWriterTest extends TestCase {

  public SmartBedWriterTest(String name) {
    super(name);
  }

  public void test() throws IOException {
    final ByteArrayOutputStream bos = new ByteArrayOutputStream();
    final SmartBedWriter c = new SmartBedWriter(bos);
    final DiscordBedRecord rec = createRecord("chr2", 123);
    final DiscordBedRecord rec1 = createRecord("chr1", 12);
    final DiscordBedRecord rec2 = createRecord("chr1", 10);
    final DiscordBedRecord rec3 = createRecord("chr1", 30);
    rec3.setFiltered();

    c.addRecord(rec);
    c.addRecord(rec1);
    c.addRecord(rec2);
    c.addRecord(rec3);

    c.close();

    final String exp = ""
      + "chr2  123 128" + LS // Note that reordering is only within chromosomes
      + "chr1  10  15" + LS
      + "chr1  12  17" + LS;

    assertEquals(exp.replaceAll("[ ]+", "\t"), TestUtils.stripLines(bos.toString(), "#", LS));

    assertTrue(bos.toString().contains("#chr1\t30\t35"));
  }

  private static DiscordBedRecord createRecord(String chr, int pos) {
    return new DiscordBedRecord(chr, pos, pos + 5);
  }

}
