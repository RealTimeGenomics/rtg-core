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

import java.util.Properties;

import com.rtg.vcf.VcfReaderTest;
import com.rtg.vcf.VcfRecord;

/**
 */
public class NullModelBuilderTest extends AbstractModelBuilderTest<NullModelBuilder> {

  @Override
  NullModelBuilder getModelBuilder(String[] formatAttributes, String[] infoAttributes, String[] derivedAttributes) {
    return new NullModelBuilder(formatAttributes, infoAttributes, derivedAttributes);
  }

  public void testAnnotate() {
    final NullModelBuilder amb = getModelBuilder(new String[]{"GQ"}, new String[]{}, new String[]{});
    assertNotNull(amb);
    amb.setModelParameters(new Properties());
    amb.build();
    final AbstractPredictModel apm = amb.getModel();
    assertNotNull(apm);
    final VcfRecord record = VcfReaderTest.vcfLineToRecord("chr5\t12041\trs55926606\tA\tT\t100\tPASS\tXRX\tGT:GQ\t0|1:12.3\t.\t1|1:99.9");
    final String s = record.toString();
    apm.annotate(record);
    assertEquals(s, record.toString());
  }

}
