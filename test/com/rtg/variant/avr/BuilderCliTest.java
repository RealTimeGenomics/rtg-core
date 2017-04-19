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
import java.io.IOException;

import com.rtg.launcher.AbstractCli;
import com.rtg.launcher.AbstractCliTest;
import com.rtg.util.TestUtils;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.TestDirectory;
import com.rtg.util.test.FileHelper;

/**
 */
public class BuilderCliTest extends AbstractCliTest {

  @Override
  protected AbstractCli getCli() {
    return new BuilderCli();
  }

  public void testHelp() {
    checkHelp("Create an AVR model",
        "--derived-annotations=STRING", "derived fields to use in model",
        "--format-annotations=STRING", "FORMAT fields to use in model",
        "--info-annotations=STRING", "INFO fields to use in model",
        "--qual-annotation", "use QUAL annotation in model",
        "--sample", "the name of the sample to select (required when using multi-sample VCF files)",
        "-n,", "--negative=FILE", "negative training examples",
        "-o,", "--output=FILE", "output AVR model",
        "-p,", "--positive=FILE", "positive training examples"
        );
  }

  public void testRun() throws IOException {
    try (final TestDirectory dir = new TestDirectory()) {
      final File posVcf = new File(dir, "pos.vcf");
      FileHelper.resourceToFile("com/rtg/variant/avr/resources/positives.vcf", posVcf);
      final File negVcf = new File(dir, "neg.vcf");
      FileHelper.resourceToFile("com/rtg/variant/avr/resources/negatives.vcf", negVcf);

      final File avr = new File(dir, "model.avr");

      checkMainInitOk("-o", avr.getPath(), "-n", negVcf.getPath(), "-p", posVcf.getPath(), "--format-annotations", "GT,GQ", "--info-annotations", "XRX,RCE", "--derived-annotations", "IC,EP");

      ModelFactory mf = new ModelFactory(avr, 0.0);
      assertTrue(mf.getModel() instanceof MlAvrPredictModel);
      assertEquals("GT,GQ", mf.getModelProperties().getProperty(AbstractModelBuilder.MODEL_PROPERTY_FORMAT_ANNOTATIONS));
      assertEquals("XRX,RCE", mf.getModelProperties().getProperty(AbstractModelBuilder.MODEL_PROPERTY_INFO_ANNOTATIONS));
      assertEquals("IC,EP", mf.getModelProperties().getProperty(AbstractModelBuilder.MODEL_PROPERTY_DERIVED_ANNOTATIONS));
      assertEquals("ML", mf.getModelProperties().getProperty(AbstractModelBuilder.MODEL_PROPERTY_TYPE));
      assertTrue(avr.delete());

      checkMainInitOk("-o", avr.getPath(), "-n", negVcf.getPath(), "-p", posVcf.getPath(), "--format-annotations", "GT,GQ", "--info-annotations", "XRX,RCE", "--derived-annotations", "IC,EP", "--sample", "NA12878", "--qual-annotation");

      mf = new ModelFactory(avr, 0.0);
      assertTrue(mf.getModel() instanceof MlAvrPredictModel);
      assertEquals("GT,GQ", mf.getModelProperties().getProperty(AbstractModelBuilder.MODEL_PROPERTY_FORMAT_ANNOTATIONS));
      assertEquals("XRX,RCE", mf.getModelProperties().getProperty(AbstractModelBuilder.MODEL_PROPERTY_INFO_ANNOTATIONS));
      assertEquals("IC,EP", mf.getModelProperties().getProperty(AbstractModelBuilder.MODEL_PROPERTY_DERIVED_ANNOTATIONS));
      assertEquals("ML", mf.getModelProperties().getProperty(AbstractModelBuilder.MODEL_PROPERTY_TYPE));
      assertTrue(mf.getModel().toString(), mf.getModel().toString().contains("DERIVED-EP(DOUBLE),DERIVED-IC(DOUBLE),FORMAT-GQ(INTEGER),FORMAT-GT(STRING),INFO-RCE(BOOLEAN),INFO-XRX(BOOLEAN),QUAL(DOUBLE)"));
      assertTrue(avr.delete());

      final File props = new File(dir, "props.txt");
      FileHelper.resourceToFile("com/rtg/variant/avr/resources/gtcomplex.props", props);

      checkMainInitOk("-o", avr.getPath(), "-n", negVcf.getPath(), "-p", posVcf.getPath(), "--qual-annotation", "--Xmodel-type", "gt-complex", "--Xmodel-params", props.getPath());

      mf = new ModelFactory(avr, 0.0);
      assertTrue(mf.getModel() instanceof GtQualComplexMultiplierModel);
      assertEquals("", mf.getModelProperties().getProperty(AbstractModelBuilder.MODEL_PROPERTY_FORMAT_ANNOTATIONS));
      assertEquals("", mf.getModelProperties().getProperty(AbstractModelBuilder.MODEL_PROPERTY_INFO_ANNOTATIONS));
      assertEquals("", mf.getModelProperties().getProperty(AbstractModelBuilder.MODEL_PROPERTY_DERIVED_ANNOTATIONS));
      assertEquals("GT_COMPLEX", mf.getModelProperties().getProperty(AbstractModelBuilder.MODEL_PROPERTY_TYPE));
      assertTrue(avr.delete());
      final GtQualComplexMultiplierModel model = (GtQualComplexMultiplierModel) mf.getModel();
      TestUtils.containsAll(model.toString(),
        "multiplier.gq.simple.homozygous\t2.1",
        "multiplier.gq.simple.heterozygous\t0.612",
        "multiplier.gq.complex.homozygous\t1.26",
        "multiplier.gq.complex.heterozygous\t0.342",
        "multiplier.qual.simple\t0.3",
        "multiplier.qual.complex\t0.5"
      );

      String error = checkMainInitBadFlags("-o", avr.getPath(), "-n", negVcf.getPath(), "-p", posVcf.getPath(), "--Xmodel-type", "gt-complex", "--Xmodel-params", props.getPath(), "--derived-annotations", "NOTDERIVED");
      TestUtils.containsAll(error, "Invalid value \"NOTDERIVED\" for flag --derived-annotations");

      error = checkMainInitBadFlags("-o", avr.getPath(), "-n", negVcf.getPath(), "-p", posVcf.getPath(), "--qual-annotation", "--Xmodel-type", "gt-complex", "--Xmodel-params", props.getPath(), "--sample", "SAMPLE");
      TestUtils.containsAll(error, "Sample name not found in VCF file: SAMPLE : " + posVcf.getPath());

      final File fakeVcf = new File(dir, "fake.vcf");
      FileUtils.stringToFile("##fileformat=VCFv4.1\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA12878\tNA12888\tNA12345\n", fakeVcf);

      error = checkMainInitBadFlags("-o", avr.getPath(), "-n", fakeVcf.getPath(), "-p", posVcf.getPath(), "--qual-annotation");
      TestUtils.containsAll(error, "Need to specify a sample name for a multi-sample VCF file: " + fakeVcf.getPath());

    }
  }

}
