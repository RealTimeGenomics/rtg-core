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

import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;

import com.rtg.vcf.VcfAnnotator;
import com.rtg.vcf.VcfRecord;
import com.rtg.vcf.header.VcfHeader;

/**
 *
 * Implementations of this abstract class should be thread safe and should not do any internal threading.
 *
 */
public abstract class AbstractPredictModel implements VcfAnnotator {

  private static final String DEFAULT_OUTPUT_FIELD_NAME = "AVR";

  private final String mField;

  /**
   * Create a new predict model from the contents of the given input stream.
   * @param is input stream to read from
   */
  public AbstractPredictModel(InputStream is) {
    this(is, DEFAULT_OUTPUT_FIELD_NAME);
  }

  /**
   * Create a new predict model from the contents of the given input stream, and write scores to the given field.
   * @param is input stream to read from.
   * @param field VCF field name to write score to.
   */
  public AbstractPredictModel(InputStream is, String field) {
    mField = field;
  }

  /**
   * Constructor used for building model from scratch with model builder mechanism.
   */
  AbstractPredictModel() {
    mField = DEFAULT_OUTPUT_FIELD_NAME;
  }

  /**
   * Process the given VCF record through the model annotating the record with a new score value.
   * @param record VCF record to process.
   */
  @Override
  public abstract void annotate(VcfRecord record);

  /**
   * Process the given VCF record through the model annotating the record for given sample with a new score value.
   * @param record VCF record to process.
   * @param sampleNumber sample to annotate.
   */
  public abstract void annotateSample(VcfRecord record, int sampleNumber);

  /**
   * Add comments to output VCF header.
   * @param header VCF header to annotate.
   */
  @Override
  public abstract void updateHeader(VcfHeader header);

  @Override
  public abstract String toString();

  /**
   * Returns a summary report about the predictions that have been made.
   * @return a report
   */
  public abstract String getSummary();

  /**
   * Saves the model to the given output stream.
   * @param os output stream to save to.
   * @throws IOException if an error occurs.
   */
  public abstract void save(OutputStream os) throws IOException;

  /**
   * Return the field name to store score in.
   * @return field name.
   */
  String getField() {
    return mField;
  }
}
