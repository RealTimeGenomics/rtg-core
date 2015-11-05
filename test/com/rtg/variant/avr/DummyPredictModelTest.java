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
import java.io.InputStreamReader;
import java.io.OutputStream;

import com.rtg.variant.avr.DummyPredictModelTest.DummyPredictModel;
import com.rtg.vcf.VcfRecord;
import com.rtg.vcf.header.MetaType;
import com.rtg.vcf.header.VcfHeader;
import com.rtg.vcf.header.VcfNumber;

/**
 */
public class DummyPredictModelTest extends AbstractPredictModelTest<DummyPredictModel> {

  protected static class DummyPredictModel extends AbstractPredictModel {
    String mData = null;

    public DummyPredictModel(InputStream is) throws IOException {
      super(is);
      try (final InputStreamReader isr = new InputStreamReader(is)) {
        final char[] buff = new char[1000];
        final int n = isr.read(buff);
        mData = new String(buff, 0, n);
      }
    }

    public DummyPredictModel() {
      mData = "123";
    }

    @Override
    public void updateHeader(VcfHeader header) {
      header.addFormatField(getField(), MetaType.FLOAT, VcfNumber.ONE, "AVR score");
    }

    @Override
    public String toString() {
      return mData;
    }

    @Override
    public void save(OutputStream os) throws IOException {
      final byte[] b = mData.getBytes();
      os.write(b);
      os.flush();
    }

    @Override
    public String getSummary() {
      return "A summary";
    }

    @Override
    public void annotate(VcfRecord record) {
      for (int s = 0; s < record.getNumberOfSamples(); s++) {
        annotateSample(record, s);
      }
    }

    @Override
    public void annotateSample(VcfRecord record, int sampleNo) {
      record.setFormatAndSample(getField(), mData, sampleNo);
    }
  }

  @Override
  DummyPredictModel getPredictModel() {
    return new DummyPredictModel();
  }

  @Override
  DummyPredictModel getPredictModel(InputStream is) throws IOException {
    return new DummyPredictModel(is);
  }
}
