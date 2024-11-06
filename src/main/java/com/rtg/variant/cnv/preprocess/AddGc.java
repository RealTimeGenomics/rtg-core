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
package com.rtg.variant.cnv.preprocess;

import java.io.IOException;
import java.util.Map;

import com.rtg.mode.DNA;
import com.rtg.reader.ReaderUtils;
import com.rtg.reader.SequencesReader;
import com.rtg.util.diagnostic.NoTalkbackSlimException;
import com.rtg.util.intervals.SequenceNameLocus;

/**
 * Computes GC and percent GC content of each region
 */
public class AddGc implements DatasetProcessor {

  static final String GC_NAME = "gc_content_abs";
  static final String PCT_GC_NAME = "gc_content_rel";

  private static final int N_BYTE = DNA.N.ordinal();
  private static final int C_BYTE = DNA.C.ordinal();
  private static final int G_BYTE = DNA.G.ordinal();

  private final SequencesReader mSequencesReader;
  private final Map<String, Long> mNames;

  /**
   * Constructor
   * @param reader reader supplying reference base information
   * @throws IOException if there is a problem reading the reference data
   */
  public AddGc(SequencesReader reader) throws IOException {
    mSequencesReader = reader;
    mNames = ReaderUtils.getSequenceNameMap(mSequencesReader);
  }

  private double computeGcPct(int gcCount, int nCount, int length) {
    final int nonN = length - nCount;
    return nonN == 0 ? 0 : (double) gcCount / nonN;
  }

  @Override
  public void process(RegionDataset dataset) throws IOException {
    final NumericColumn gcCol = dataset.addColumn(new IntColumn(GC_NAME));
    final NumericColumn pgcCol = dataset.addColumn(new NumericColumn(PCT_GC_NAME));
    String chr = null;
    final byte[] template = new byte[(int) mSequencesReader.maxLength()];
    int length = 0;
    for (int j = 0; j < dataset.regions().size(); ++j) {
      final SequenceNameLocus r = dataset.regions().get(j);
      if (!r.getSequenceName().equals(chr)) {
        chr = r.getSequenceName();
        if (!mNames.containsKey(chr)) {
          throw new NoTalkbackSlimException("Reference SDF does not contain sequence '" + chr + "'");
        }
        length = mSequencesReader.read(mNames.get(chr), template);
      }
      if (r.getEnd() > length) {
        throw new NoTalkbackSlimException("Region is out of range for the supplied chromosome! " + r);
      }
      int gc = 0;
      int nc = 0;
      for (int i = r.getStart(); i < r.getEnd(); ++i) {
        final byte base = template[i];
        if (base == C_BYTE || base == G_BYTE) {
          ++gc;
        } else if (base == N_BYTE) {
          ++nc;
        }
      }
      gcCol.add(gc);
      pgcCol.add(computeGcPct(gc, nc, r.getLength()));
    }
  }
}
