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
import java.util.Properties;
import java.util.zip.ZipEntry;
import java.util.zip.ZipInputStream;

import com.rtg.util.diagnostic.NoTalkbackSlimException;


/**
 * Load AVR models from an AVR zip file.
 *
 */
public class ModelFactory {

  final double mAvrThreshold;
  AbstractPredictModel mModel = null;
  final Properties mProperties = new Properties();

  /**
   * Loads a model from the given AVR model file.  The correct model class is determined from the contents of the file.
   *
   * @param avrFile AVR model file
   * @throws IOException if an error occurs reading the file
   */
  public ModelFactory(File avrFile) throws IOException {
    this(avrFile, 0);
  }

  /**
   * Loads a model from the given AVR model file.  The correct model class is determined from the contents of the file.
   *
   * @param avrFile AVR model file
   * @param avrThreshold fail records below this threshold
   * @throws IOException if an error occurs reading the file
   */
  public ModelFactory(File avrFile, double avrThreshold) throws IOException {
    mAvrThreshold = avrThreshold;

    // load appropriate AVRPredictModel based on the contents of the file
    // 1) load properties from model.properties in zip file
    // 2) from header info determine appropriate predict model
    // 3) create predict model with model file content

    try (final ZipInputStream zin = new ZipInputStream(new FileInputStream(avrFile))) {
      ZipEntry ze = zin.getNextEntry();
      if (ze == null || !AbstractModelBuilder.MODEL_PROPERTIES_FILE_NAME.equals(ze.getName())) {
        throw new NoTalkbackSlimException("FIle does not look like an AVR model: " + avrFile.getPath());
      }

      //System.err.println(">>> " + ze.getName());
      mProperties.load(zin);

      // check version number from properties
      final String version = mProperties.getProperty(AbstractModelBuilder.MODEL_AVR_VERSION);
      if (version == null) {
        throw new NoTalkbackSlimException("No version number in AVR model: " + avrFile.getPath());
      } else {
        final int avrVersion = Integer.parseInt(version);
        if (avrVersion != AbstractModelBuilder.AVR_VERSION) {
          throw new NoTalkbackSlimException("Cannot handle AVR model with a version of " + avrVersion);
        }
      }

      ze = zin.getNextEntry();
      if (ze == null || !AbstractModelBuilder.MODEL_FILE_NAME.equals(ze.getName())) {
        throw new NoTalkbackSlimException("FIle does not look like an AVR model: " + avrFile.getPath());
      }

      // get model class to load from properties
      final String modelTypeName = mProperties.getProperty(AbstractModelBuilder.MODEL_PROPERTY_TYPE);
      if (modelTypeName == null) {
        // no model name in AVR header
        throw new NoTalkbackSlimException("No model type specied in AVR file.");
      }
      final ModelType modelType;
      try {
        modelType = ModelType.valueOf(modelTypeName);
      } catch (IllegalArgumentException iae) {
        // unknown model name in AVR header
        throw new NoTalkbackSlimException("Unknown AVR model type: " + modelTypeName);
      }
      // add other models to load in this switch statement
      try {
        switch (modelType) {
          case GT_COMPLEX:
            mModel = new GtQualComplexMultiplierModel(zin);
            break;
          case ML:
            mModel = new MlAvrPredictModel(zin, mAvrThreshold);
            break;
          case NULL:
            mModel = new NullModel(zin);
            break;
          default:
            throw new UnsupportedOperationException("Loading support not implemented for AVR model type: " + modelType);
        }
      } catch (RuntimeException re) {
        throw new NoTalkbackSlimException(re, "Could not load AVR model: " + avrFile.getPath() + ". It may be created with a newer version of RTG, or may be corrupt.");
      }
    }
  }


  public AbstractPredictModel getModel() {
    return mModel;
  }

  public Properties getModelProperties() {
    return mProperties;
  }
}
