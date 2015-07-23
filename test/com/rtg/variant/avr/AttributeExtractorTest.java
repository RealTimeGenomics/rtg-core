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

import static com.rtg.util.StringUtils.LS;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import com.rtg.ml.Dataset;
import com.rtg.ml.Instance;
import com.rtg.ml.MlDataType;
import com.rtg.util.TestUtils;
import com.rtg.util.io.TestDirectory;
import com.rtg.util.test.FileHelper;
import com.rtg.variant.avr.AttributeExtractor.IncompatibleHeaderException;
import com.rtg.vcf.VcfReader;
import com.rtg.vcf.annotation.AnnotationDataType;
import com.rtg.vcf.header.MetaType;

import junit.framework.TestCase;

/**
 */
public class AttributeExtractorTest extends TestCase {
  public void testConstructor() {
    try {
      new AttributeExtractor((Annotation) null);
      fail("accepted null annotations");
    } catch (NullPointerException npe) {
      // expected
    }

    try {
      new AttributeExtractor(new FormatAnnotation("ABC", AnnotationDataType.STRING), null);
      fail("accepted null annotations");
    } catch (NullPointerException npe) {
      // expected
    }
  }

  public void testCheckHeader() throws IOException {
    try (final TestDirectory dir = new TestDirectory()) {
      final File vcf = new File(dir, "a.vcf");
      FileHelper.resourceToFile("com/rtg/variant/avr/resources/positives.vcf", vcf);

      try (final VcfReader reader = VcfReader.openVcfReader(vcf)) {
        final List<Annotation> annotations = new ArrayList<>();
        annotations.add(new QualAnnotation());
        annotations.add(new InfoAnnotation("XRX", AnnotationDataType.BOOLEAN));
        annotations.add(new InfoAnnotation("RCE", AnnotationDataType.BOOLEAN));
        annotations.add(new InfoAnnotation("CT", AnnotationDataType.INTEGER));
        annotations.add(new FormatAnnotation("GT", AnnotationDataType.STRING));
        annotations.add(new FormatAnnotation("AB", AnnotationDataType.DOUBLE));
        annotations.add(new FormatAnnotation("DP", AnnotationDataType.INTEGER));
        AttributeExtractor ae = new AttributeExtractor(annotations.toArray(new Annotation[annotations.size()]));
        try {
          ae.checkHeader(reader.getHeader());
        } catch (IncompatibleHeaderException ihe) {
          fail("Got exception with valid header");
        }

        annotations.add(0, new InfoAnnotation("ABC", AnnotationDataType.BOOLEAN)); // make first annotation in list...
        annotations.add(new InfoAnnotation("BAD1", AnnotationDataType.INTEGER));
        annotations.add(new InfoAnnotation("BAD2", AnnotationDataType.INTEGER));
        annotations.add(new InfoAnnotation("BAD4", AnnotationDataType.DOUBLE));
        annotations.add(new FormatAnnotation("GQ", AnnotationDataType.STRING));
        annotations.add(new FormatAnnotation("DEF", AnnotationDataType.INTEGER));
        annotations.add(new FormatAnnotation("RS", AnnotationDataType.STRING));
        annotations.add(new FormatAnnotation("BAD3", AnnotationDataType.DOUBLE));
        ae = new AttributeExtractor(annotations.toArray(new Annotation[annotations.size()]));
        try {
          ae.checkHeader(reader.getHeader());
          fail("Failed to get exception with invalid header");
        } catch (IncompatibleHeaderException ihe) {
          final String message = ihe.getMessage();
          //System.err.println(message);
          TestUtils.containsAll(message,
              "INFO field does not exist: ABC",
              "INFO field type mismatch: BAD4 : Integer != DOUBLE",
              "INFO field arbitary number of values: BAD1",
              "INFO field too many values: BAD2 : 2 > 1",
              "FORMAT field type mismatch: GQ : Integer != STRING",
              "FORMAT field does not exist: DEF",
              "FORMAT field arbitary number of values: RS",
              "FORMAT field too many values: BAD3 : 3 != 1"
              );
        }
      }
    }
  }

  public void testGetDataset() throws IOException {
    try (final TestDirectory dir = new TestDirectory()) {
      final File vcf = new File(dir, "a.vcf");
      FileHelper.resourceToFile("com/rtg/variant/avr/resources/positives.vcf", vcf);

      try (final VcfReader reader = VcfReader.openVcfReader(vcf)) {
        final List<Annotation> annotations = new ArrayList<>();
        annotations.add(new QualAnnotation());
        annotations.add(new InfoAnnotation("XRX", AnnotationDataType.BOOLEAN));
        annotations.add(new InfoAnnotation("RCE", AnnotationDataType.BOOLEAN));
        annotations.add(new InfoAnnotation("CT", AnnotationDataType.INTEGER));
        annotations.add(new FormatAnnotation("GT", AnnotationDataType.STRING));
        annotations.add(new FormatAnnotation("AB", AnnotationDataType.DOUBLE));
        annotations.add(new FormatAnnotation("DP", AnnotationDataType.INTEGER));
        final AttributeExtractor ae = new AttributeExtractor(annotations.toArray(new Annotation[annotations.size()]));
        try {
          ae.checkHeader(reader.getHeader());
        } catch (IncompatibleHeaderException ihe) {
          System.err.println(ihe.getMessage());
          fail("Got exception with valid header");
        }

        final Dataset ds = ae.getDataset();
        assertNotNull(ds);

        reader.hasNext();
        final double[] instance = ae.getInstance(reader.next(), 0);
        assertNotNull(instance);
        assertEquals(7, instance.length);
//        assertTrue(instance[0] instanceof Double); // AB
//        assertTrue(instance[1] instanceof Integer); // DP
//        assertTrue(instance[2] instanceof String);  // GT
        assertTrue(Double.isNaN(instance[3])); // CT
//        assertTrue(instance[4] instanceof Boolean); // RCE
//        assertTrue(instance[5] instanceof Boolean); // XRX
//        assertTrue(instance[6] instanceof Double); // QUAL
        assertEquals(0.564, (double) ds.getAttributes()[0].decodeValue(instance[0]), 0.01);
        assertEquals(39, (int) ds.getAttributes()[1].decodeValue(instance[1]));
        assertEquals("1/0", (String) ds.getAttributes()[2].decodeValue(instance[2]));
        assertFalse((boolean) ds.getAttributes()[4].decodeValue(instance[4]));
        assertFalse((boolean) ds.getAttributes()[5].decodeValue(instance[5]));
        assertEquals(50.3, (double) ds.getAttributes()[6].decodeValue(instance[6]), 0.01);

        ds.addInstance(new Instance(instance, true));
      }
    }
  }

  public void testMultisampleGetDataset() throws IOException {
    try (final TestDirectory dir = new TestDirectory()) {
      final File vcf = new File(dir, "a.vcf");
      FileHelper.resourceToFile("com/rtg/variant/avr/resources/multisample.vcf", vcf);

      try (final VcfReader reader = VcfReader.openVcfReader(vcf)) {
        final List<Annotation> annotations = new ArrayList<>();
        annotations.add(new QualAnnotation());
        annotations.add(new InfoAnnotation("XRX", AnnotationDataType.BOOLEAN));
        annotations.add(new InfoAnnotation("RCE", AnnotationDataType.BOOLEAN));
        annotations.add(new InfoAnnotation("DP", AnnotationDataType.INTEGER));
        annotations.add(new FormatAnnotation("GT", AnnotationDataType.STRING));
        annotations.add(new FormatAnnotation("AB", AnnotationDataType.DOUBLE));
        annotations.add(new FormatAnnotation("DP", AnnotationDataType.INTEGER));
        final AttributeExtractor ae = new AttributeExtractor(annotations.toArray(new Annotation[annotations.size()]));
        try {
          ae.checkHeader(reader.getHeader());
        } catch (IncompatibleHeaderException ihe) {
          System.err.println(ihe.getMessage());
          fail("Got exception with valid header");
        }

        final Dataset ds = ae.getDataset();
        assertNotNull(ds);

        reader.hasNext();
        double[] instance = ae.getInstance(reader.next(), 1);
        assertNotNull(instance);
        assertEquals(7, instance.length);
//        assertTrue(instance[0] instanceof Double); // AB
//        assertTrue(instance[1] instanceof Integer); // DP
//        assertTrue(instance[2] instanceof String);  // GT
//        assertTrue(instance[3] instanceof Integer); // DP info
//        assertTrue(instance[4] instanceof Boolean); // RCE
//        assertTrue(instance[5] instanceof Boolean); // XRX
//        assertTrue(instance[6] instanceof Double); // QUAL

        assertEquals(0.444, (double) ds.getAttributes()[0].decodeValue(instance[0]), 0.01);
        assertEquals(18, (int) ds.getAttributes()[1].decodeValue(instance[1]));
        assertEquals("1/0", (String) ds.getAttributes()[2].decodeValue(instance[2]));
        assertEquals(66, (int) ds.getAttributes()[3].decodeValue(instance[3]));
        assertFalse((boolean) ds.getAttributes()[4].decodeValue(instance[4]));
        assertFalse((boolean) ds.getAttributes()[5].decodeValue(instance[5]));
        assertEquals(319.2, (double) ds.getAttributes()[6].decodeValue(instance[6]), 0.01);

        ds.addInstance(new Instance(instance, true));

        reader.hasNext();
        instance = ae.getInstance(reader.next(), 2);
        assertNotNull(instance);
        assertEquals(7, instance.length);
        assertTrue(Double.isNaN(instance[0])); // AB
//        assertTrue(instance[1] instanceof Integer); // DP
//        assertTrue(instance[2] instanceof String);  // GT
//        assertTrue(instance[3] instanceof Integer); // DP info
//        assertTrue(instance[4] instanceof Boolean); // RCE
//        assertTrue(instance[5] instanceof Boolean); // XRX
//        assertTrue(instance[6] instanceof Double); // QUAL

        assertEquals(6, (int) ds.getAttributes()[1].decodeValue(instance[1]));
        assertEquals("1/1", (String) ds.getAttributes()[2].decodeValue(instance[2]));
        assertEquals(18, (int) ds.getAttributes()[3].decodeValue(instance[3]));
        assertFalse((boolean) ds.getAttributes()[4].decodeValue(instance[4]));
        assertTrue((boolean) ds.getAttributes()[5].decodeValue(instance[5]));
        assertEquals(1823.4, (double) ds.getAttributes()[6].decodeValue(instance[6]), 0.01);

        ds.addInstance(new Instance(instance, false));

        assertEquals("Number of examples with missing values:" + LS
            + "  FORMAT-AB  1" + LS
            + "  FORMAT-DP  0" + LS
            + "  FORMAT-GT  0" + LS
            + "  INFO-DP    0" + LS
            + "  INFO-RCE   0" + LS
            + "  INFO-XRX   0" + LS
            + "  QUAL       0" + LS
            , ae.missingValuesReport());
      }
    }
  }

  public void testGetCompatibleType() {
    assertEquals(AnnotationDataType.INTEGER, AttributeExtractor.getCompatibleType(MetaType.INTEGER));
    assertEquals(AnnotationDataType.STRING, AttributeExtractor.getCompatibleType(MetaType.CHARACTER));
    assertEquals(AnnotationDataType.BOOLEAN, AttributeExtractor.getCompatibleType(MetaType.FLAG));
    assertEquals(AnnotationDataType.DOUBLE, AttributeExtractor.getCompatibleType(MetaType.FLOAT));
    assertEquals(AnnotationDataType.STRING, AttributeExtractor.getCompatibleType(MetaType.STRING));
  }

  public void testGetArffAttributeText() {
    assertEquals("@attribute FORMAT-XXX numeric", AttributeExtractor.getArffAttributeText(new FormatAnnotation("XXX", AnnotationDataType.INTEGER)));
    assertEquals("@attribute INFO-YYY {" + Boolean.TRUE + "," + Boolean.FALSE + "}", AttributeExtractor.getArffAttributeText(new InfoAnnotation("YYY", AnnotationDataType.BOOLEAN)));
    assertEquals("@attribute FORMAT-XXX numeric", AttributeExtractor.getArffAttributeText(new FormatAnnotation("XXX", AnnotationDataType.DOUBLE)));
    assertEquals("@attribute FORMAT-XXX string", AttributeExtractor.getArffAttributeText(new FormatAnnotation("XXX", AnnotationDataType.STRING)));
  }

  public void testGetValueAsArffString() {
    assertEquals("?", AttributeExtractor.getValueAsArffString(null));
    assertEquals("2.1", AttributeExtractor.getValueAsArffString(2.1));
    assertEquals("abc", AttributeExtractor.getValueAsArffString("abc"));
  }

  public void testGetMlDataType() {
    assertEquals(MlDataType.STRING, AttributeExtractor.getMlDataType(AnnotationDataType.STRING));
    assertEquals(MlDataType.DOUBLE, AttributeExtractor.getMlDataType(AnnotationDataType.DOUBLE));
    assertEquals(MlDataType.INTEGER, AttributeExtractor.getMlDataType(AnnotationDataType.INTEGER));
    assertEquals(MlDataType.BOOLEAN, AttributeExtractor.getMlDataType(AnnotationDataType.BOOLEAN));
  }
}
