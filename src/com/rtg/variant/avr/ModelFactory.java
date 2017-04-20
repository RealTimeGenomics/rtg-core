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
