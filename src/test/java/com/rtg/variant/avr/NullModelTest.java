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
package com.rtg.variant.avr;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.InputStream;

import com.rtg.vcf.VcfReaderTest;
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

    final VcfRecord record = VcfReaderTest.vcfLineToRecord("chr5\t12041\trs55926606\tA\tT\t100\tPASS\tXRX\tGT:GQ\t0|1:12.3\t.\t1|1:99.9");
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
