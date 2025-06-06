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

import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;

import com.rtg.launcher.globals.GlobalFlags;
import com.rtg.launcher.globals.CoreGlobalFlags;
import com.rtg.ml.MlPredictLoader;
import com.rtg.ml.PredictClassifier;
import com.rtg.util.StringUtils;
import com.rtg.util.Utils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.NoTalkbackSlimException;
import com.rtg.vcf.VcfRecord;
import com.rtg.vcf.header.FilterField;
import com.rtg.vcf.header.FormatField;
import com.rtg.vcf.header.MetaType;
import com.rtg.vcf.header.VcfHeader;
import com.rtg.vcf.header.VcfNumber;

/**
 * Machine learning predict model.
 */
public class MlAvrPredictModel extends AbstractPredictModel {

  static final int SERIAL_VERSION = 1;

  private AttributeExtractor mAttributeExtractor = null;
  private PredictClassifier mClassifier = null;

  private static final int NUM_BINS = 20;
  private final long[] mScoreBins = new long[NUM_BINS + 1];
  int mCurrentVersion;
  private String mFilterName;
  final double mPredictionThreshold;


  /**
   * Create a new predict model from the contents of the given input stream.
   * @param is input stream to read from
   * @param threshold min AVR score threshold
   * @throws IOException if error occurs loading model
   */
  public MlAvrPredictModel(InputStream is, double threshold) throws IOException {
    load(is);
    mPredictionThreshold = threshold;
    assert mClassifier != null;
    assert mAttributeExtractor != null;
  }

  private void load(InputStream is) throws IOException {
    // load attribute extractor
    // then load classifier
    final DataInputStream dis = new DataInputStream(is);
    mCurrentVersion = dis.readInt();
    if (mCurrentVersion == 1) {
      mAttributeExtractor = AttributeExtractor.load(is);
      mClassifier = MlPredictLoader.loadPredictClassifier(is, mAttributeExtractor.getDataset());
    } else {
      throw new IOException("Unsupported model version: " + mCurrentVersion);
    }
  }

  /**
   * Constructor used for building model from scratch with model builder mechanism.
   * @param classifier the classifier to wrap
   */
  MlAvrPredictModel(PredictClassifier classifier) {
    super();
    if (classifier == null) {
      throw new NullPointerException();
    }
    mClassifier = classifier;
    mPredictionThreshold = 0;
  }

  void setAttributeExtractor(AttributeExtractor ae) {
    mAttributeExtractor = ae;
  }

  @Override
  public void annotate(VcfRecord record) {
    boolean aboveThreshold = false;
    for (int s = 0; s < record.getNumberOfSamples(); ++s) {
      final double prediction = annotateSampleNoPadding(record, s);
      if (prediction >= mPredictionThreshold) {
        aboveThreshold = true;
      }
    }
    if (!aboveThreshold) {
      record.addFilter(mFilterName);
    }
  }

  @Override
  public void annotateSample(VcfRecord record, int sampleNumber) {
    annotateSampleNoPadding(record, sampleNumber);
    record.padFormatAndSample(getField());
  }

  private double annotateSampleNoPadding(VcfRecord record, int sampleNumber) {
    // extract fields from record to build instance object array
    final double[] instance = mAttributeExtractor.getInstance(record, sampleNumber);
    final double prediction = mClassifier.predict(instance);
    record.setFormatAndSample(getField(), Utils.realFormat(prediction, 4), sampleNumber);
    incrementScore(prediction);
    return prediction;
  }


  @Override
  public void updateHeader(VcfHeader header) {
    // Check compatibility with the attributes that are expected
    try {
      mAttributeExtractor.checkHeader(header);
    } catch (AttributeExtractor.IncompatibleHeaderException ihe) {
      if (GlobalFlags.isSet(CoreGlobalFlags.AVR_ALLOW_UNDECLARED_ATTRIBUTES)) {
        Diagnostic.warning("VCF does not declare all fields in AVR model, they will be treated as missing values:" + StringUtils.LS + ihe.getMessage());
      } else {
        throw new NoTalkbackSlimException("The input VCF header is missing required fields:" + StringUtils.LS + ihe.getMessage());
      }
    }
    mFilterName = "AVR" + mPredictionThreshold;
    if (mPredictionThreshold > 0) {
      header.ensureContains(new FilterField(mFilterName, "AVR score below " + mPredictionThreshold));
    }
    header.ensureContains(new FormatField(getField(), MetaType.FLOAT, VcfNumber.ONE, "AVR score"));
  }

  @Override
  public String toString() {
    final StringBuilder sb = new StringBuilder();
    sb.append("ATTRIBUTE EXTRACTOR:").append(StringUtils.LS);
    if (mAttributeExtractor == null) {
      sb.append("NOT SET").append(StringUtils.LS);
    } else {
      sb.append(mAttributeExtractor).append(StringUtils.LS);
    }
    sb.append("MODEL:").append(StringUtils.LS);
    mClassifier.toString(sb, "", mAttributeExtractor != null ? mAttributeExtractor.getDataset() : null);
    return sb.toString();
  }

  private void incrementScore(double value) {
    mScoreBins[(int) (NUM_BINS * value)]++;
  }

  @Override
  public String getSummary() {
    final StringBuilder sb = new StringBuilder();
    for (int i = 0; i < NUM_BINS; ++i) {
      sb.append("[")
      .append(Utils.realFormat((double) i / NUM_BINS, 2))
      .append("..")
      .append(Utils.realFormat((double) (i + 1) / NUM_BINS, 2))
      .append(")\t")
      .append(mScoreBins[i])
      .append(StringUtils.LS);
    }
    sb.append("1\t")
    .append(mScoreBins[NUM_BINS])
    .append(StringUtils.LS);
    return sb.toString();
  }

  @Override
  public void save(OutputStream os) throws IOException {
    final DataOutputStream dos = new DataOutputStream(os);
    dos.writeInt(SERIAL_VERSION);
    // save attribute extractor
    mAttributeExtractor.save(os);
    // save classifier
    mClassifier.save(dos, mAttributeExtractor.getDataset());
  }

}
