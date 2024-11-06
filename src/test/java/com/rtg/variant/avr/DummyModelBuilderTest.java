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
