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
