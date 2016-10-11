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
    return new Properties(mProperties);
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
    try (final ZipOutputStream zout = new ZipOutputStream(new FileOutputStream(avrFile))) {
      zout.setLevel(9); // maximum compression

      // save properties
      zout.putNextEntry(new ZipEntry(MODEL_PROPERTIES_FILE_NAME));
      mProperties.store(zout, MODEL_PROPERTIES_COMMENT);

      // save model
      zout.putNextEntry(new ZipEntry(MODEL_FILE_NAME));
      getModel().save(zout);
    }
  }

}
