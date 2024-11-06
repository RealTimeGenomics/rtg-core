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
