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
import java.io.FileOutputStream;
import java.io.IOException;
import java.text.SimpleDateFormat;
import java.util.Arrays;
import java.util.Date;
import java.util.Properties;
import java.util.UUID;
import java.util.zip.ZipEntry;
import java.util.zip.ZipOutputStream;

import com.rtg.util.StringUtils;
import com.rtg.util.cli.CommandLine;


/**
 * Builder associated with an {@link AbstractPredictModel}.
 *
 */
public abstract class AbstractModelBuilder<T extends AbstractPredictModel> {
  /** AVR version */
  public static final String MODEL_AVR_VERSION = "avr.version";
  /** AVR Model type */
  public static final String MODEL_PROPERTY_TYPE = "type";
  /** Date model was created */
  public static final String MODEL_PROPERTY_DATE = "date";
  /** Unique model identifier */
  public static final String MODEL_PROPERTY_MODEL_ID = "id";
  /** List of FORMAT fields used in model */
  public static final String MODEL_PROPERTY_FORMAT_ANNOTATIONS = "format.annotations";
  /** List of INFO fields used in model */
  public static final String MODEL_PROPERTY_INFO_ANNOTATIONS = "info.annotations";
  /** List of derived fields used in model */
  public static final String MODEL_PROPERTY_DERIVED_ANNOTATIONS = "derived.annotations";
  /** Whether QUAL field is used in model */
  public static final String MODEL_PROPERTY_QUAL_ANNOTATION = "qual.annotation.used";

  /** Command line used to generate model */
  public static final String MODEL_PROPERTY_COMMAND_LINE = "command.line";

  private static final String MODEL_PROPERTIES_COMMENT = "RTG AVR Model";
  /** Internal model properties file name */
  public static final String MODEL_PROPERTIES_FILE_NAME = "model.properties";
  /** Internal model file name */
  public static final String MODEL_FILE_NAME = "model";


  /** Current version */
  public static final int AVR_VERSION = 1;

  final Properties mProperties = new Properties();
  final Properties mParameters = new Properties();

  Boolean mUseQualAttribute = Boolean.FALSE;
  final String[] mFormatAttributes;
  final String[] mInfoAttributes;
  final String[] mDerivedAttributes;

  T mModel = null;

  /**
   * Create an AVR model by extracting values from the given attribute fields.
   * @param formatAttributes FORMAT field attributes to use
   * @param infoAttributes INFO field attributes to use
   * @param derivedAttributes derived field attributes to use
   */
  public AbstractModelBuilder(String[] formatAttributes, String[] infoAttributes, String[] derivedAttributes) {
    for (String x : formatAttributes) {
      if (x == null) {
        throw new NullPointerException();
      }
    }
    for (String x : infoAttributes) {
      if (x == null) {
        throw new NullPointerException();
      }
    }
    for (String x : derivedAttributes) {
      if (x == null) {
        throw new NullPointerException();
      }
    }
    mFormatAttributes = Arrays.copyOf(formatAttributes, formatAttributes.length);
    mInfoAttributes = Arrays.copyOf(infoAttributes, infoAttributes.length);
    mDerivedAttributes = Arrays.copyOf(derivedAttributes, derivedAttributes.length);

    mProperties.setProperty(MODEL_AVR_VERSION, Integer.toString(AVR_VERSION));
    mProperties.setProperty(MODEL_PROPERTY_MODEL_ID, UUID.randomUUID().toString());
    mProperties.setProperty(MODEL_PROPERTY_DATE, new SimpleDateFormat("yyyy-MM-dd-HH-mm-ss").format(new Date()));
    mProperties.setProperty(MODEL_PROPERTY_FORMAT_ANNOTATIONS, StringUtils.join(",", mFormatAttributes));
    mProperties.setProperty(MODEL_PROPERTY_INFO_ANNOTATIONS, StringUtils.join(",", mInfoAttributes));
    mProperties.setProperty(MODEL_PROPERTY_DERIVED_ANNOTATIONS, StringUtils.join(",", mDerivedAttributes));
    mProperties.setProperty(MODEL_PROPERTY_QUAL_ANNOTATION, Boolean.toString(mUseQualAttribute));
    String cmd = CommandLine.getCommandLine();
    if (cmd == null) {
      cmd = "NO COMMAND LINE";
    }
    mProperties.setProperty(MODEL_PROPERTY_COMMAND_LINE, cmd);

  }

  /**
   * Set model specific parameters. Must be called before build to ensure parameters are used.
   * @param parameters parameters and properties.
   */
  public void setModelParameters(Properties parameters) {
    // model specific properties - read from --XX flags?
    mParameters.putAll(parameters);
  }

  /**
   * Set whether the AVR module should use the QUAL attribute.
   * @param flag whether to use QUAL
   */
  public void useQualAttribute(Boolean flag) {
    mUseQualAttribute = flag;
    mProperties.setProperty(MODEL_PROPERTY_QUAL_ANNOTATION, Boolean.toString(mUseQualAttribute));
  }

  /**
   * Build a module using the given VCF datasets as training data.
   * @param vcfDatasets an array of datasets
   * @throws IOException if an error occurs reading the files
   */
  public abstract void build(VcfDataset... vcfDatasets) throws IOException;

  /**
   * Return model properties for inclusion in the AVR file.
   *
   * @return model properties properties
   */
  public Properties getModelPropeties() {
    // check properties have be set???
    if (mProperties.getProperty(MODEL_PROPERTY_TYPE) == null) {
      throw new IllegalStateException();
    }
    return mProperties;
  }

  /**
   * Return the model that has been built.
   * @return AVR prediction model.
   */
  public T getModel() {
    return mModel;
  }

  /**
   * Writes the properties and model to the given file.
   *
   * @param avrFile file to write to
   * @throws IOException if a problem occurs while writing the file
   */
  public void save(File avrFile) throws IOException {
    save(avrFile, mProperties, getModel());
  }

  static <U extends AbstractPredictModel> void save(File avrFile, Properties properties, U model) throws IOException {
    try (final ZipOutputStream zout = new ZipOutputStream(new FileOutputStream(avrFile))) {
      zout.setLevel(9); // maximum compression

      // save properties
      zout.putNextEntry(new ZipEntry(MODEL_PROPERTIES_FILE_NAME));
      properties.store(zout, MODEL_PROPERTIES_COMMENT);

      // save model
      zout.putNextEntry(new ZipEntry(MODEL_FILE_NAME));
      model.save(zout);
    }
  }
}
