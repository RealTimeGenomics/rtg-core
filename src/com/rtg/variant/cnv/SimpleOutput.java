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
package com.rtg.variant.cnv;


import java.io.IOException;
import java.io.OutputStream;

import com.rtg.util.StringUtils;
import com.rtg.util.Utils;
import com.rtg.util.integrity.Exam;
import com.rtg.util.integrity.IntegralAbstract;

/**
 */
public class SimpleOutput extends IntegralAbstract implements LSMOutput {

  private static final byte[] CNV = "cnv".getBytes();
  private static final byte[] LSB = StringUtils.LS.getBytes();

  private final OutputStream mOut;
  private final int mBucketSize;

  /**
   * @param out output stream
   * @param bucketSize bucket size
   */
  public SimpleOutput(final OutputStream out, int bucketSize) {
    mOut = out;
    mBucketSize = bucketSize;
  }

  @Override
  public void out(final String id, final int start, final int end, final double sum, final double sum2) throws IOException {
    final int length = end - start;
    final double mean = sum / length;
    final double threshold = sum * mean;
    final double stdDev = length == 1 ? 0.0 : Math.sqrt((sum2 - threshold) / (length - 1));
    mOut.write(id.getBytes());
    mOut.write('\t');
    mOut.write(Integer.toString(ntPos(start)).getBytes());
    mOut.write('\t');
    mOut.write(Integer.toString(ntPos(end)).getBytes());
    mOut.write('\t');
    mOut.write(CNV);
    mOut.write('\t');
    mOut.write(Utils.realFormat(mean, 3).getBytes());
    mOut.write('\t');
    mOut.write(Utils.realFormat(stdDev, 3).getBytes());
    mOut.write(LSB);
  }

  private int ntPos(int blockPos) {
    return (blockPos - 1) * mBucketSize;
  }

  @Override
  public void close() throws IOException {
    mOut.close();
  }

  @Override
  public void toString(final StringBuilder sb) {
    sb.append("SimpleOutput");
  }

  @Override
  public boolean integrity() {
    Exam.assertTrue(mOut != null);
    return true;
  }

}
