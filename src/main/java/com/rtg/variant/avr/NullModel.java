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
import java.io.InputStream;
import java.io.OutputStream;
import java.util.Properties;

import com.rtg.vcf.VcfRecord;
import com.rtg.vcf.header.VcfHeader;

/**
 * A model that does nothing.
 */
public class NullModel extends AbstractPredictModel {

  private static final String PROPERTIES_COMMENT = "Null model.";

  /**
   * Default constructor.
   */
  public NullModel() {
  }

  /**
   * Creates a model from the model file.
   *
   * @param is input stream to read from
   * @throws IOException if an error occurs loading from stream
   */
  public NullModel(final InputStream is) throws IOException {
    super();
    load(is);
  }

  private void load(InputStream is) throws IOException {
    // Just load the properties even though they are ignored
    new Properties().load(is);
  }

  @Override
  public void save(OutputStream os) throws IOException {
    new Properties().store(os, PROPERTIES_COMMENT);
  }

  @Override
  public void annotate(VcfRecord record) {
    // do nothing
  }

  @Override
  public void annotateSample(VcfRecord record, int sampleNo) {
    // do nothing
  }

  @Override
  public void updateHeader(VcfHeader header) {
    // Nothing to add
  }

  @Override
  public String toString() {
    return "";
  }

  @Override
  public String getSummary() {
    return "";
  }

}
