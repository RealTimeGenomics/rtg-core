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

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.InputStream;

import com.rtg.vcf.VcfReader;
import com.rtg.vcf.VcfRecord;
import com.rtg.vcf.header.VcfHeader;

/**
 */
public class NullModelTest extends AbstractPredictModelTest<NullModel> {

  @Override
  NullModel getPredictModel() {
    return new NullModel();
  }

  @Override
  NullModel getPredictModel(InputStream is) throws IOException {
    return new NullModel(is);
  }

  @Override
  public void testLoadSave() throws Exception {
    final AbstractPredictModel apm = getPredictModel();
    assertNotNull(apm);
    assertEquals("AVR", apm.getField());
    assertNotNull(apm.toString());
    try (final ByteArrayOutputStream baos = new ByteArrayOutputStream()) {
      apm.save(baos);
      try (final ByteArrayInputStream bais = new ByteArrayInputStream(baos.toByteArray())) {
        final AbstractPredictModel apm2 = getPredictModel(bais);
        assertNotNull(apm2);
        assertEquals("AVR", apm2.getField());
        assertEquals(apm.toString(), apm2.toString());
      }
    }
  }

  @Override
  public void testAnnotate() {
    final AbstractPredictModel apm = getPredictModel();
    assertNotNull(apm);
    assertEquals("AVR", apm.getField());

    final VcfRecord record = VcfReader.vcfLineToRecord("chr5\t12041\trs55926606\tA\tT\t100\tPASS\tXRX\tGT:GQ\t0|1:12.3\t.\t1|1:99.9");
    final String s = record.toString();
    apm.annotate(record);
    assertEquals(s, record.toString());
    apm.annotate(record);
    assertEquals(s, record.toString());
  }

  @Override
  public void testUpdateHeader() {
    final VcfHeader header = new VcfHeader();
    final AbstractPredictModel apm = getPredictModel();
    assertNotNull(apm);
    apm.updateHeader(header);
    assertEquals(0, header.getInfoLines().size());
    assertEquals(0, header.getFormatLines().size());
  }

}
