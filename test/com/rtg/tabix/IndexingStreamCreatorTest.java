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
package com.rtg.tabix;

import java.io.File;
import java.io.IOException;

import com.rtg.util.io.TestDirectory;
import com.rtg.vcf.VcfRecord;
import com.rtg.vcf.VcfWriter;
import com.rtg.vcf.header.VcfHeader;

import junit.framework.TestCase;

/**
 */
public class IndexingStreamCreatorTest extends TestCase {

  public void testIt() throws IOException {
    try (TestDirectory dir = new TestDirectory("isct")) {
      final File vcfFile = new File(dir, "thingy.vcf.gz");
      VcfHeader header = new VcfHeader();
      header.addCommonHeader();
      header.addSampleName("sample");
      try (IndexingStreamCreator streamHandler = new IndexingStreamCreator(vcfFile, System.out, true, new TabixIndexer.VcfIndexerFactory(), true)) {
        try (VcfWriter writer = new VcfWriter(header, streamHandler.createStreamsAndStartThreads())) {
          VcfRecord rec = new VcfRecord();
          rec.setNumberOfSamples(1);
          rec.setSequence("foo");
          rec.setStart(1000);
          rec.setRefCall("A");
          rec.addAltCall("T");
          rec.addFormatAndSample("GT", "1/0");
          writer.write(rec);
        }
      }
      assertTrue(new File(dir, "thingy.vcf.gz.tbi").exists());
    }
  }
}
