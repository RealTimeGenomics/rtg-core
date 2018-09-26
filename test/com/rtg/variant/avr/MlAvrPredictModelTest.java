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

import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.util.Properties;

import com.rtg.ml.Attribute;
import com.rtg.ml.BuildClassifier;
import com.rtg.ml.BuilderFactory;
import com.rtg.ml.Dataset;
import com.rtg.ml.Instance;
import com.rtg.ml.MlDataType;
import com.rtg.ml.RandomTreeBuilder;
import com.rtg.ml.ZeroRBuilder;
import com.rtg.util.Resources;
import com.rtg.util.diagnostic.NoTalkbackSlimException;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.TestDirectory;
import com.rtg.util.test.FileHelper;
import com.rtg.vcf.VcfReader;
import com.rtg.vcf.VcfRecord;
import com.rtg.vcf.header.VcfHeader;

/**
 */
public class MlAvrPredictModelTest extends AbstractPredictModelTest<MlAvrPredictModel> {
  @Override
  MlAvrPredictModel getPredictModel() {
    final Dataset ds = new Dataset(new Attribute("empty", MlDataType.DOUBLE));
    ds.addInstance(new Instance(new double[1], true));
    ds.addInstance(new Instance(new double[1], false));
    final BuildClassifier builder = new ZeroRBuilder();
    builder.build(ds);
    final MlAvrPredictModel model = new MlAvrPredictModel(builder.getClassifier());
    model.setAttributeExtractor(new AttributeExtractor(new QualAnnotation()));
    return model;
  }

  @Override
  MlAvrPredictModel getPredictModel(InputStream is) throws IOException {
    return new MlAvrPredictModel(is, 0.0);
  }

  public void testAnnotate2() throws IOException {
    final String[] formatAttributes = {"GQ"};
    final String[] infoAttributes = {"XRX"};
    final String[] derivedAttributes = {};
    final MlAvrModelBuilder amb = new MlAvrModelBuilder(formatAttributes, infoAttributes, derivedAttributes);
    try (final TestDirectory dir = new TestDirectory()) {
      final File posVcf = new File(dir, "pos.vcf");
      FileHelper.resourceToFile("com/rtg/variant/avr/resources/positives.vcf", posVcf);
      final File negVcf = new File(dir, "neg.vcf");
      FileHelper.resourceToFile("com/rtg/variant/avr/resources/negatives.vcf", negVcf);

      final Properties params = new Properties();
      //params.setProperty(MlAvrModelBuilder.PARAMETER_ML_MODEL_TYPE, BuilderFactory.BuilderType.BAGGED.name());
      //params.setProperty(MlAvrModelBuilder.PARAMETER_ML_MODEL_TYPE, BuilderFactory.BuilderType.ZERO_R.name());
      params.setProperty(MlAvrModelBuilder.PARAMETER_ML_MODEL_TYPE, BuilderFactory.BuilderType.RANDOM_TREE.name());
      params.setProperty(RandomTreeBuilder.PROP_MIN_INSTANCES, "" + 1);
      amb.setModelParameters(params);

      amb.build(
          new VcfDataset(posVcf, 0, VcfDataset.Classifications.ALL_POSITIVE, true, 1.0),
          new VcfDataset(negVcf, 0, VcfDataset.Classifications.ALL_NEGATIVE, true, 1.0)
      );

      final File model = new File(dir, "model.avr");
      amb.save(model);

      final ModelFactory fact = new ModelFactory(model, 0.0);
      final AbstractPredictModel apm = fact.getModel();

      assertNotNull(apm);
      assertEquals("AVR", apm.getField());

      try {
        apm.updateHeader(new VcfHeader());
        fail();
      } catch (NoTalkbackSlimException e) {
        // Expected
      }

      // Annotate all samples in the record
      VcfRecord record = VcfReader.vcfLineToRecord("chr5\t12041\trs55926606\tA\tT\t100\tPASS\tXRX\tGT:GQ\t0|1:12\t.\t1|1:99\t0|0:34");
      apm.annotate(record);

      //System.err.println(record.toString());
      assertTrue(record.toString().contains(apm.getField()));
      assertEquals(1.0, Double.valueOf(record.getFormat(apm.getField()).get(0)), 0.001);
      assertEquals(0.5127, Double.valueOf(record.getFormat(apm.getField()).get(1)), 0.001);
      assertEquals(1.0, Double.valueOf(record.getFormat(apm.getField()).get(2)), 0.001);
      assertEquals(0.0, Double.valueOf(record.getFormat(apm.getField()).get(3)), 0.001);


      // Annotate just requested samples
      record = VcfReader.vcfLineToRecord("chr5\t12041\trs55926606\tA\tT\t100\tPASS\tXRX\tGT:GQ\t0|1:12\t.\t1|1:99\t0|0:34");
      apm.annotateSample(record, 3);
      apm.annotateSample(record, 2);

      //System.err.println(record.toString());
      assertTrue(record.toString().contains(apm.getField()));
      assertEquals(".", record.getFormat(apm.getField()).get(0));
      assertEquals(".", record.getFormat(apm.getField()).get(1));
      assertEquals(1.0, Double.valueOf(record.getFormat(apm.getField()).get(2)), 0.001);
      assertEquals(0.0, Double.valueOf(record.getFormat(apm.getField()).get(3)), 0.001);
    }
  }

  private static final int MIN_VERSION = 1;
  public void testLoadVersionX() throws IOException {
    for (int i = MIN_VERSION; i <= MlAvrPredictModel.SERIAL_VERSION; ++i) {
      checkLoadVersion(i);
    }
  }

  private void checkLoadVersion(int version) throws IOException {
    final InputStream is = Resources.getResourceAsStream("com/rtg/variant/avr/resources/testMlAvrPredictModelVersion_" + version);
    try (final DataInputStream dis = new DataInputStream(is)) {
      final MlAvrPredictModel bs = new MlAvrPredictModel(dis, 0.0);
      assertEquals(version, bs.mCurrentVersion);
      final String str = bs.toString();
      assertTrue(str.contains("44/99") && str.contains("QUAL(DOUBLE)") && str.contains("0R:"));
    }
  }

  private static MlAvrPredictModel createTestModel() {
    final MlAvrPredictModel predictModel = new MlAvrPredictModel(new ZeroRBuilder.ZeroRClassifier(44, 55));
    predictModel.setAttributeExtractor(new AttributeExtractor(new QualAnnotation()));
    return predictModel;
  }

  /**
   * Creates test serialized form for current serial version
   * @param args directory in which test data should be saved
   * @throws IOException it happens
   */
  public static void main(String[] args) throws IOException {
    final File dir = new File(args[0]);
    final File output = new File(dir, "testMlAvrPredictModelVersion_" + MlAvrPredictModel.SERIAL_VERSION);
    try (DataOutputStream dos = new DataOutputStream(FileUtils.createOutputStream(output))) {
      createTestModel().save(dos);
    }
  }
}
