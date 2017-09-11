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

  /** Constant used to denote AVR scores in VCF output. */
  public static final String AVR = "AVR";

  private String mField = AVR;

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

  /**
   * Set the name of the field where the score is stored
   * @param field the field name
   */
  void setField(String field) {
    mField = field;
  }
}
