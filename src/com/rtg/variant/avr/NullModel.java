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
    super(is);
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
