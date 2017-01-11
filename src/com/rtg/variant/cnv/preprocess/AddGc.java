/*
 * Copyright (c) 2016. Real Time Genomics Limited.
 *
 * Use of this source code is bound by the Real Time Genomics Limited Software Licence Agreement
 * for Academic Non-commercial Research Purposes only.
 *
 * If you did not receive a license accompanying this file, a copy must first be obtained by email
 * from support@realtimegenomics.com.  On downloading, using and/or continuing to use this source
 * code you accept the terms of that license agreement and any amendments to those terms that may
 * be made from time to time by Real Time Genomics Limited.
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
   * @param reader supplying reference base information
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
        throw new NoTalkbackSlimException("Region is out of range for the supplied chromosome! " + r.toString());
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
