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
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.Properties;
import java.util.zip.ZipEntry;
import java.util.zip.ZipOutputStream;

import com.rtg.util.TestUtils;
import com.rtg.util.cli.CommandLine;
import com.rtg.util.diagnostic.NoTalkbackSlimException;
import com.rtg.util.io.TestDirectory;
import com.rtg.util.test.FileHelper;

import junit.framework.TestCase;

/**
 */
public class ModelFactoryTest extends TestCase {

  public void testConstructor() throws IOException {
    CommandLine.clearCommandArgs();
    try {
      new ModelFactory(null, 0.0);
      fail("accepted null file");
    } catch (NullPointerException npe) {
      // expected
    }

    try (TestDirectory dir = new TestDirectory()) {
      final File file = new File(dir, "model.avr");
      FileHelper.stringToGzFile("hello", file);
      try {
        new ModelFactory(file, 0.0);
        fail("accepted bad file");
      } catch (NoTalkbackSlimException ntse) {
        //System.err.println(ntse.getMessage()); // expected
        assertTrue(ntse.getMessage().contains("FIle does not look like an AVR model"));
      }
    }

    try {
      new ModelFactory(null, 0.0);
      fail("accepted null file");
    } catch (NullPointerException npe) {
      // expected
    }

    try (TestDirectory dir = new TestDirectory()) {
      final File file = new File(dir, "model.avr");
      try (final ZipOutputStream zout = new ZipOutputStream(new FileOutputStream(file))) {
        zout.setLevel(9); // maximum compression
        // save properties
        zout.putNextEntry(new ZipEntry("XXXX"));
      }
      try {
        new ModelFactory(file, 0.0);
        fail("accepted bad file");
      } catch (NoTalkbackSlimException ntse) {
        //System.err.println(ntse.getMessage()); // expected
        assertTrue(ntse.getMessage().contains("FIle does not look like an AVR model"));
      }
    }


    try (TestDirectory dir = new TestDirectory()) {
      final File file = new File(dir, "model.avr");
      try (final ZipOutputStream zout = new ZipOutputStream(new FileOutputStream(file))) {
        zout.setLevel(9); // maximum compression
        // save properties
        zout.putNextEntry(new ZipEntry(AbstractModelBuilder.MODEL_PROPERTIES_FILE_NAME));
      }
      try {
        new ModelFactory(file, 0.0);
        fail("accepted bad file");
      } catch (NoTalkbackSlimException ntse) {
        //System.err.println(ntse.getMessage()); // expected
        assertTrue(ntse.getMessage().contains("No version number in AVR model"));
      }
    }

    try (TestDirectory dir = new TestDirectory()) {
      final File file = new File(dir, "model.avr");
      try (final ZipOutputStream zout = new ZipOutputStream(new FileOutputStream(file))) {
        zout.setLevel(9); // maximum compression
        // save properties
        zout.putNextEntry(new ZipEntry(AbstractModelBuilder.MODEL_PROPERTIES_FILE_NAME));
        final Properties props = new Properties();
        props.put(AbstractModelBuilder.MODEL_AVR_VERSION, "2");
        props.store(zout, "?");
      }
      try {
        new ModelFactory(file, 0.0);
        fail("accepted bad file");
      } catch (NoTalkbackSlimException ntse) {
        //System.err.println(ntse.getMessage()); // expected
        assertTrue(ntse.getMessage().contains("Cannot handle AVR model with a version of 2"));
      }
    }

    try (TestDirectory dir = new TestDirectory()) {
      final File file = new File(dir, "model.avr");
      try (final ZipOutputStream zout = new ZipOutputStream(new FileOutputStream(file))) {
        zout.setLevel(9); // maximum compression
        // save properties
        zout.putNextEntry(new ZipEntry(AbstractModelBuilder.MODEL_PROPERTIES_FILE_NAME));
        final Properties props = new Properties();
        props.put(AbstractModelBuilder.MODEL_AVR_VERSION, "1");
        props.store(zout, "?");

        zout.putNextEntry(new ZipEntry("XXXX"));
      }
      try {
        new ModelFactory(file, 0.0);
        fail("accepted bad file");
      } catch (NoTalkbackSlimException ntse) {
        //System.err.println(ntse.getMessage()); // expected
        assertTrue(ntse.getMessage().contains("FIle does not look like an AVR model"));
      }
    }

    final String[] formatAttributes = {"GP", "DP", "RE", "AB"};
    final String[] infoAttributes = {"XRX", "RCE", "SP"};
    final String[] derivedAttributes = {"IC", "EP"};
    final AbstractModelBuilder<?> amb = new GtQualComplexMultiplierModelBuilder(formatAttributes, infoAttributes, derivedAttributes);
    assertNotNull(amb);
    assertNull(amb.getModel());

    amb.build();

    try (TestDirectory dir = new TestDirectory()) {
      final File file = new File(dir, "model.avr");
      amb.save(file);

      final ModelFactory amf = new ModelFactory(file, 0.0);

      final AbstractPredictModel apm = amf.getModel();
      assertNotNull(apm);
      assertTrue(apm instanceof GtQualComplexMultiplierModel);
      assertEquals("AVR", apm.getField());
      TestUtils.containsAll(apm.toString(),
          "GQ multipliers:",
          "multiplier.gq.simple.homozygous\t2.0",
          "multiplier.gq.simple.heterozygous\t0.512",
          "multiplier.gq.complex.homozygous\t1.16",
          "multiplier.gq.complex.heterozygous\t0.242",
          "QUAL multipliers:",
          "multiplier.qual.simple\t0.2",
          "multiplier.qual.complex\t0.04"
          );

      final Properties props = amf.getModelProperties();
      assertEquals("1", props.getProperty(AbstractModelBuilder.MODEL_AVR_VERSION));
      assertEquals("GP,DP,RE,AB", props.getProperty(AbstractModelBuilder.MODEL_PROPERTY_FORMAT_ANNOTATIONS));
      assertEquals("XRX,RCE,SP", props.getProperty(AbstractModelBuilder.MODEL_PROPERTY_INFO_ANNOTATIONS));
      assertEquals("IC,EP", props.getProperty(AbstractModelBuilder.MODEL_PROPERTY_DERIVED_ANNOTATIONS));
      assertNotNull(props.getProperty(AbstractModelBuilder.MODEL_PROPERTY_DATE));
      assertNotNull(props.getProperty(AbstractModelBuilder.MODEL_PROPERTY_MODEL_ID));
      assertEquals("NO COMMAND LINE", props.getProperty(AbstractModelBuilder.MODEL_PROPERTY_COMMAND_LINE));

    }
  }

  public void testModelLoading() throws IOException {
    try (TestDirectory dir = new TestDirectory()) {
      final File file = new File(dir, "model.avr");
      FileHelper.resourceToFile("com/rtg/variant/avr/resources/gtqualcomplexmodel.avr", file);

      ModelFactory amf = new ModelFactory(file, 0.0);

      AbstractPredictModel apm = amf.getModel();
      assertNotNull(apm);
      assertTrue(apm instanceof GtQualComplexMultiplierModel);
      assertTrue(file.delete());
      Properties props = amf.getModelProperties();
      assertEquals("GT_COMPLEX", props.getProperty(AbstractModelBuilder.MODEL_PROPERTY_TYPE));
      assertEquals("1", props.getProperty(AbstractModelBuilder.MODEL_AVR_VERSION));

      FileHelper.resourceToFile("com/rtg/variant/avr/resources/mlmodel.avr", file);

      amf = new ModelFactory(file, 0.0);

      apm = amf.getModel();
      assertNotNull(apm);
      assertTrue(apm instanceof MlAvrPredictModel);
      assertTrue(file.delete());
      props = amf.getModelProperties();
      assertEquals("ML", props.getProperty(AbstractModelBuilder.MODEL_PROPERTY_TYPE));
      assertEquals("1", props.getProperty(AbstractModelBuilder.MODEL_AVR_VERSION));
    }
  }

  public void testNullModelLoading() throws IOException {
    try (TestDirectory dir = new TestDirectory()) {
      final File file = new File(dir, "model.avr");
      FileHelper.resourceToFile("com/rtg/variant/avr/resources/null.avr", file);
      final ModelFactory amf = new ModelFactory(file, 0.0);
      final AbstractPredictModel apm = amf.getModel();
      assertNotNull(apm);
      assertTrue(apm instanceof NullModel);
      assertTrue(file.delete());
    }
  }
}
