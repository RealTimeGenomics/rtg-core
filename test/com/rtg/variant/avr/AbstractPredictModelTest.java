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

package com.rtg.variant.avr;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.List;

import com.rtg.util.io.TestDirectory;
import com.rtg.vcf.VcfReader;
import com.rtg.vcf.VcfRecord;
import com.rtg.vcf.header.MetaType;
import com.rtg.vcf.header.VcfHeader;

import junit.framework.TestCase;

/**
 * Abstract unit tests for {@link AbstractPredictModel}.
 */
public abstract class AbstractPredictModelTest<T extends AbstractPredictModel> extends TestCase {
  abstract T getPredictModel();
  abstract T getPredictModel(InputStream is) throws IOException;

  public void testLoadSave() throws Exception {
    final AbstractPredictModel apm = getPredictModel();
    assertNotNull(apm);
    assertEquals("AVR", apm.getField());
    assertNotNull(apm.toString());

    try (TestDirectory dir = new TestDirectory()) {
      final File file = File.createTempFile("model", "asr", dir);
      try (final FileOutputStream fos = new FileOutputStream(file)) {
        apm.save(fos);
      }

      try (final FileInputStream fis = new FileInputStream(file)) {
        final AbstractPredictModel apm2 = getPredictModel(fis);
        assertNotNull(apm2);
        assertEquals("AVR", apm2.getField());
        assertEquals(apm.toString(), apm2.toString());
      }
    }
  }

  public void testAnnotate() {
    final AbstractPredictModel apm = getPredictModel();
    assertNotNull(apm);
    assertEquals("AVR", apm.getField());

    final VcfRecord record = VcfReader.vcfLineToRecord("chr5\t12041\trs55926606\tA\tT\t100\tPASS\tXRX\tGT:GQ\t0|1:12.3\t.\t1|1:99.9");
    apm.annotate(record);

    //System.err.println(record.toString());
    assertTrue(record.toString().contains("AVR"));
    final List<String> avrVals = record.getFormatAndSample().get("AVR");
    assertNotNull(avrVals);
    assertEquals(3, avrVals.size());
  }

  public void testUpdateHeader() {
    final VcfHeader header = new VcfHeader();
    final AbstractPredictModel apm = getPredictModel();
    assertNotNull(apm);
    assertEquals("AVR", apm.getField());
    apm.updateHeader(header);
    //System.err.println(header);
    assertEquals(1, header.getFormatLines().size());
    assertEquals("AVR", header.getFormatLines().get(0).getId());
    assertEquals(1, header.getFormatLines().get(0).getNumber().getNumber());
    assertEquals(MetaType.FLOAT, header.getFormatLines().get(0).getType());
  }
}
