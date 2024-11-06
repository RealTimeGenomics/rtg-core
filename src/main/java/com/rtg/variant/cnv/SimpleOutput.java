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
package com.rtg.variant.cnv;


import java.io.IOException;
import java.io.OutputStream;

import com.rtg.util.StringUtils;
import com.rtg.util.Utils;
import com.rtg.util.integrity.Exam;
import com.rtg.util.integrity.IntegralAbstract;

/**
 */
public class SimpleOutput extends IntegralAbstract implements LsmOutput {

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
    final double stdDev = length == 1 ? 0.0 : Math.sqrt((sum2 - sum * mean) / (length - 1));
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
