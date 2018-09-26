/*
 * Copyright (c) 2017. Real Time Genomics Limited.
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

import com.rtg.ml.Attribute;
import com.rtg.ml.MlDataType;
import com.rtg.util.io.TestDirectory;
import com.rtg.util.test.FileHelper;

/**
 */
public class ArffDatasetWriterTest extends AbstractModelBuilderTest<ArffDatasetWriter> {

  @Override
  ArffDatasetWriter getModelBuilder(String[] formatAttributes, String[] infoAttributes, String[] derivedAttributes) {
    return new ArffDatasetWriter(formatAttributes, infoAttributes, derivedAttributes);
  }

  public void testGetArffAttributeText() {
    assertEquals("@attribute FORMAT-XXX numeric", ArffDatasetWriter.getArffAttributeText(new FormatAnnotation("XXX", AnnotationDataType.INTEGER), new Attribute("XXX", MlDataType.INTEGER)));
    assertEquals("@attribute INFO-YYY {" + Boolean.TRUE + "," + Boolean.FALSE + "}", ArffDatasetWriter.getArffAttributeText(new InfoAnnotation("YYY", AnnotationDataType.BOOLEAN), new Attribute("XXX", MlDataType.BOOLEAN)));
    assertEquals("@attribute FORMAT-XXX numeric", ArffDatasetWriter.getArffAttributeText(new FormatAnnotation("XXX", AnnotationDataType.DOUBLE), new Attribute("XXX", MlDataType.DOUBLE)));
    final Attribute sAtt = new Attribute("XXX", MlDataType.STRING);
    sAtt.encodeValue("a1");
    sAtt.encodeValue("a2");
    assertEquals("@attribute FORMAT-XXX {a1,a2}", ArffDatasetWriter.getArffAttributeText(new FormatAnnotation("XXX", AnnotationDataType.STRING), sAtt));
    for (int i = 3; i < ArffDatasetWriter.MAX_NOMINAL + 4; i++) {
      sAtt.encodeValue("a" + i);
    }
    assertEquals("@attribute FORMAT-XXX string", ArffDatasetWriter.getArffAttributeText(new FormatAnnotation("XXX", AnnotationDataType.STRING), sAtt));
  }

  public void testGetValueAsArffString() {
    assertEquals("?", ArffDatasetWriter.getValueAsArffString(null));
    assertEquals("2.1", ArffDatasetWriter.getValueAsArffString(2.1));
    assertEquals("abc", ArffDatasetWriter.getValueAsArffString("abc"));
  }

  @Override
  public void testLoadSave() throws Exception {
    final String[] formatAttributes = {"GQ", "DP", "RE", "AB"};
    final String[] infoAttributes = {"XRX", "RCE"};
    final String[] derivedAttributes = {"IC", "EP", "RA", "QD"};
    final ArffDatasetWriter amb = getModelBuilder(formatAttributes, infoAttributes, derivedAttributes);
    assertNotNull(amb);
    assertNull(amb.getModel());

    try (final TestDirectory dir = new TestDirectory()) {
      final File posVcf = FileHelper.resourceToFile("com/rtg/variant/avr/resources/positives.vcf", new File(dir, "pos.vcf"));
      final File negVcf = FileHelper.resourceToFile("com/rtg/variant/avr/resources/negatives.vcf", new File(dir, "neg.vcf"));

      amb.build(
        new VcfDataset(posVcf, 0, VcfDataset.Classifications.ALL_POSITIVE, false, 2.0),
        new VcfDataset(negVcf, 0, VcfDataset.Classifications.ALL_NEGATIVE, false, 2.0)
      );

      final File file = new File(dir, "model.arff");
      amb.save(file);
      mNano.check("arff-dataset.arff", FileHelper.fileToString(file));

      final File bothVcf = FileHelper.resourceToFile("com/rtg/variant/avr/resources/posandneg.vcf", new File(dir, "posandneg.vcf"));
      amb.build(new VcfDataset(bothVcf, 0, VcfDataset.Classifications.ANNOTATED, false, 2.0));
      final File bothfile = new File(dir, "modelboth.arff");
      amb.save(bothfile);
      mNano.check("arff-dataset.arff", FileHelper.fileToString(bothfile).replaceAll("modelboth", "model"));
    }
  }
}
