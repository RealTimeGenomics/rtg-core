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

import java.util.Properties;

import com.rtg.vcf.VcfReader;
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
    final VcfRecord record = VcfReader.vcfLineToRecord("chr5\t12041\trs55926606\tA\tT\t100\tPASS\tXRX\tGT:GQ\t0|1:12.3\t.\t1|1:99.9");
    final String s = record.toString();
    apm.annotate(record);
    assertEquals(s, record.toString());
  }

}
