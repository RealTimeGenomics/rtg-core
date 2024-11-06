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

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.List;

import com.rtg.util.io.TestDirectory;
import com.rtg.vcf.VcfReaderTest;
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

    final VcfRecord record = VcfReaderTest.vcfLineToRecord("chr5\t12041\trs55926606\tA\tT\t100\tPASS\tXRX\tGT:GQ\t0|1:12.3\t.\t1|1:99.9");
    apm.annotate(record);

    //System.err.println(record.toString());
    assertTrue(record.toString().contains("AVR"));
    final List<String> avrVals = record.getFormat("AVR");
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
