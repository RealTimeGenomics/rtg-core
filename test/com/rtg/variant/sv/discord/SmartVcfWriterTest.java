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
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.io.MemoryPrintStream;
import com.rtg.vcf.VcfRecord;
import com.rtg.vcf.header.VcfHeader;

import junit.framework.TestCase;


/**
 */
public class SmartVcfWriterTest extends TestCase {

  public SmartVcfWriterTest(String name) {
    super(name);
  }

  public void test() throws IOException {
    final ByteArrayOutputStream bos = new ByteArrayOutputStream();
    final VcfHeader header = new VcfHeader();
    header.addSampleName("SAMPLE");
    final SmartVcfWriter c = new SmartVcfWriter(bos, header);
    final VcfRecord rec = createRecord("chr2", 123);
    final VcfRecord rec1 = createRecord("chr1", 12);
    final VcfRecord rec2 = createRecord("chr1", 10).addAltCall("g"); //should be after rec3 since it has more alts
    final VcfRecord rec3 = createRecord("chr1", 10);


    c.addRecord(rec);
    c.addRecord(rec1);
    c.addRecord(rec2);
    c.addRecord(rec3);

    c.close();

    final String exp = ""
      + "chr2  123 . a c,t 12.8  . ." + LS // Note that reordering is only within chromosomes
      + "chr1  10  . a c,t 12.8  . ." + LS
      + "chr1  10  . a c,t,g 12.8  . ." + LS
      + "chr1  12  . a c,t 12.8  . ." + LS;

    assertEquals(exp.replaceAll("[ ]+", "\t"), TestUtils.stripLines(bos.toString(), "#", LS));
  }

  private static VcfRecord createRecord(String chr, int pos) {
    final VcfRecord rec = new VcfRecord();
    rec.setSequence(chr)
    .setStart(pos - 1)
    .setId(".")
    .setQuality("12.8")
    .setRefCall("a")
    .addAltCall("c")
    .addAltCall("t");
    return rec;
  }

  public void testLimitExceeded() throws IOException {
    final MemoryPrintStream mps = orderLimit(10000);
    TestUtils.containsAll(mps.toString()
        , "VcfRecord dropped due to excessive out-of-order processing"
        , "chr2\t1\t.\ta\tc,t\t12.8\t.\t."
    );
  }

  private MemoryPrintStream orderLimit(int number) throws IOException {
    final MemoryPrintStream mps = new MemoryPrintStream();
    Diagnostic.setLogStream(mps.printStream());
    try {
      final ByteArrayOutputStream bos = new ByteArrayOutputStream();
      final VcfHeader header = new VcfHeader();
      header.addSampleName("SAMPLE");
      final SmartVcfWriter c = new SmartVcfWriter(bos, header);
      // Why +3 ?
      c.addRecord(createRecord("chr2", 2)); //first output.
      c.addRecord(createRecord("chr2", number + 3)); //if number >  buffer length will cause first record to be written
      c.addRecord(createRecord("chr2", 1)); // before first output, should cause error if number is big enough
      c.close();

    } finally {
      Diagnostic.setLogStream();
    }
    return mps;
  }

  public void testUnderLimit() throws IOException {
    final MemoryPrintStream mps = orderLimit(9999);
    assertEquals("", mps.toString());
  }
}
