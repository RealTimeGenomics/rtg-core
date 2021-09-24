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
import java.io.FileInputStream;
import java.io.IOException;
import java.util.zip.ZipEntry;
import java.util.zip.ZipInputStream;

import com.rtg.variant.avr.DummyModelBuilderTest.DummyModelBuilder;
import com.rtg.variant.avr.DummyPredictModelTest.DummyPredictModel;

/**
 * Unit tests for ModelBuilder.
 *
 */
public class DummyModelBuilderTest extends AbstractModelBuilderTest<DummyModelBuilder> {

  static class DummyModelBuilder extends AbstractModelBuilder<DummyPredictModel> {

    DummyModelBuilder(String[] formatAttributes, String[] infoAttributes, String[] derivedAttributes) {
      super(formatAttributes, infoAttributes, derivedAttributes);
      mProperties.setProperty(AbstractModelBuilder.MODEL_PROPERTY_TYPE, "NULL");
    }

    @Override
    public void build(VcfDataset... vcfDatasets) {
      mModel = new DummyPredictModel();
    }

  }

  @Override
  DummyModelBuilder getModelBuilder(String[] formatAttributes, String[] infoAttributes, String[] derivedAttributes) {
    return new DummyModelBuilder(formatAttributes, infoAttributes, derivedAttributes);
  }

  static class DummyModelFactory extends ModelFactory {
    DummyModelFactory(File avrFile) throws IOException {
      super(avrFile, 0.0);
      try (final ZipInputStream zin = new ZipInputStream(new FileInputStream(avrFile))) {
        ZipEntry ze = zin.getNextEntry();
        assert ze != null;
        assert AbstractModelBuilder.MODEL_PROPERTIES_FILE_NAME.equals(ze.getName());
        //System.err.println(">>> " + ze.getName());
        mProperties.load(zin);

        //System.err.println("Properties: " + mProperties);

        assert ModelType.NULL.equals(ModelType.valueOf(mProperties.getProperty(AbstractModelBuilder.MODEL_PROPERTY_TYPE)));

        ze = zin.getNextEntry();
        assert ze != null;
        assert AbstractModelBuilder.MODEL_FILE_NAME.equals(ze.getName());
        //System.err.println(">>> " + ze.getName());
        mModel = new DummyPredictModel(zin);
      }
    }

  }

  @Override
  ModelFactory getModelFactory(File file) throws IOException {
    return new DummyModelFactory(file);
  }

}
