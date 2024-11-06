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
import java.io.IOException;
import java.util.Properties;

import com.rtg.launcher.AbstractNanoTest;
import com.rtg.util.io.TestDirectory;
import com.rtg.util.test.FileHelper;

/**
 *
 */
public abstract class AbstractModelBuilderTest<T extends AbstractModelBuilder<?>> extends AbstractNanoTest {
  abstract T getModelBuilder(String[] formatAttributes, String[] infoAttributes, String[] derivedAttributes);

  ModelFactory getModelFactory(File file) throws IOException {
    // need to override this for dummy tests
    return new ModelFactory(file, 0.0);
  }

  public void testBadConstructor() {
    final String[] nullAtts = {null};
    final String[] okAtts = {"ABC"};
    try {
      getModelBuilder(nullAtts, okAtts, okAtts);
      fail("accepted null format attribute");
    } catch (NullPointerException npe) {
      // expected
    }
    try {
      getModelBuilder(okAtts, nullAtts, okAtts);
      fail("accepted null info attribute");
    } catch (NullPointerException npe) {
      // expected
    }
    try {
      getModelBuilder(okAtts, okAtts, nullAtts);
      fail("accepted null derived attribute");
    } catch (NullPointerException npe) {
      // expected
    }
  }

  public void testConstructor() {
    final String[] formatAttributes = {"GP", "DP", "RE", "AB", "DP"};
    final String[] infoAttributes = {"XRX", "RCE", "SP", "RCE"};
    final String[] derivedAttributes = {"IC", "EP"};
    final T amb = getModelBuilder(formatAttributes, infoAttributes, derivedAttributes);
    assertNotNull(amb);

    Properties props = amb.getModelPropeties();
    assertEquals("1", props.getProperty(AbstractModelBuilder.MODEL_AVR_VERSION));
    assertEquals("GP,DP,RE,AB,DP", props.getProperty(AbstractModelBuilder.MODEL_PROPERTY_FORMAT_ANNOTATIONS));
    assertEquals("XRX,RCE,SP,RCE", props.getProperty(AbstractModelBuilder.MODEL_PROPERTY_INFO_ANNOTATIONS));
    assertEquals("IC,EP", props.getProperty(AbstractModelBuilder.MODEL_PROPERTY_DERIVED_ANNOTATIONS));
    assertEquals(Boolean.FALSE.toString(), props.getProperty(AbstractModelBuilder.MODEL_PROPERTY_QUAL_ANNOTATION));
    assertNotNull(props.getProperty(AbstractModelBuilder.MODEL_PROPERTY_DATE));
    assertNotNull(props.getProperty(AbstractModelBuilder.MODEL_PROPERTY_MODEL_ID));
    assertEquals("NO COMMAND LINE", props.getProperty(AbstractModelBuilder.MODEL_PROPERTY_COMMAND_LINE));

    amb.useQualAttribute(true);

    props = amb.getModelPropeties();
    assertEquals("1", props.getProperty(AbstractModelBuilder.MODEL_AVR_VERSION));
    assertEquals("GP,DP,RE,AB,DP", props.getProperty(AbstractModelBuilder.MODEL_PROPERTY_FORMAT_ANNOTATIONS));
    assertEquals("XRX,RCE,SP,RCE", props.getProperty(AbstractModelBuilder.MODEL_PROPERTY_INFO_ANNOTATIONS));
    assertEquals("IC,EP", props.getProperty(AbstractModelBuilder.MODEL_PROPERTY_DERIVED_ANNOTATIONS));
    assertEquals(Boolean.TRUE.toString(), props.getProperty(AbstractModelBuilder.MODEL_PROPERTY_QUAL_ANNOTATION));
    assertNotNull(props.getProperty(AbstractModelBuilder.MODEL_PROPERTY_DATE));
    assertNotNull(props.getProperty(AbstractModelBuilder.MODEL_PROPERTY_MODEL_ID));
    assertEquals("NO COMMAND LINE", props.getProperty(AbstractModelBuilder.MODEL_PROPERTY_COMMAND_LINE));


    assertNull(amb.getModel());

    try {
      ModelType.valueOf(props.getProperty(AbstractModelBuilder.MODEL_PROPERTY_TYPE));
    } catch (Exception e) {
      fail("could not encode model type to enum entry: " + AbstractModelBuilder.MODEL_PROPERTY_TYPE);
    }
  }

  public void testLoadSave() throws Exception {
    final String[] formatAttributes = {"GQ", "DP", "RE", "AB"};
    final String[] infoAttributes = {"XRX", "RCE"};
    final String[] derivedAttributes = {"IC", "EP"};
    final T amb = getModelBuilder(formatAttributes, infoAttributes, derivedAttributes);
    assertNotNull(amb);
    assertNull(amb.getModel());

    try (final TestDirectory dir = new TestDirectory()) {
      final File posVcf = FileHelper.resourceToFile("com/rtg/variant/avr/resources/positives.vcf", new File(dir, "pos.vcf"));
      final File negVcf = FileHelper.resourceToFile("com/rtg/variant/avr/resources/negatives.vcf", new File(dir, "neg.vcf"));

      amb.build(
          new VcfDataset(posVcf, 0, VcfDataset.Classifications.ALL_POSITIVE, false, 1.0),
          new VcfDataset(negVcf, 0, VcfDataset.Classifications.ALL_NEGATIVE, false, 1.0)
      );

      final File file = new File(dir, "model.avr");
      amb.save(file);

      final AbstractPredictModel apm = amb.getModel();
      checkModel(apm, getModelFactory(file).getModel());

      final File bothVcf = FileHelper.resourceToFile("com/rtg/variant/avr/resources/posandneg.vcf", new File(dir, "posandneg.vcf"));
      amb.build(new VcfDataset(bothVcf, 0, VcfDataset.Classifications.ANNOTATED, false, 1.0));

      final File bothfile = new File(dir, "modelboth.avr");
      amb.save(bothfile);
      checkModel(apm, getModelFactory(bothfile).getModel());
    }
  }

  private void checkModel(AbstractPredictModel apm, AbstractPredictModel model) {
    assertNotNull(model);
    assertEquals("AVR", model.getField());
    assertEquals(apm.toString(), model.toString());
    assertEquals(apm.getClass(), model.getClass());
  }
}
