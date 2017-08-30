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
package com.rtg.variant.cnv.region;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.Reader;
import java.util.HashMap;
import java.util.Map;
import java.util.SortedSet;
import java.util.TreeSet;

import com.rtg.mode.DNA;
import com.rtg.mode.SequenceType;
import com.rtg.reader.SequencesReader;
import com.rtg.util.array.ArrayUtils;
import com.rtg.variant.cnv.LsmOutput;
import com.rtg.variant.cnv.LeastSquaresModel;

/**
 */
public final class RegionUtils {

  private static final int MAX_UNKNOWN = DNA.N.ordinal();

  /**
   */
  enum State {
    IN, OUT
  }

  private RegionUtils() { } //prevent instantiation

  /**
   * Create a region object by reading from a reader with the format for each line:<br>
   *  name start end<br>
   *  The name is ignored and the start and end are included in the regions.
   * @param reader from the input is taken from.
   * @return the region object (never null).
   * @throws IOException whenever.
   */
  public static Map<String, Region> regionsFromReader(final Reader reader) throws IOException {
    final BufferedReader in = new BufferedReader(reader);
    final Map<String, SortedSet<AbstractCnvRegion>> map = new HashMap<>();
    while (true) {
      final String line = in.readLine();
      if (line == null) {
        break;
      }
      final String[] split = line.split("\\s+");
      assert split.length == 3 : line;
      final String seqId = split[0];
      final String sStart = split[1];
      final String sEnd = split[2];
      final SortedSet<AbstractCnvRegion> st;
      if (map.containsKey(seqId)) {
        st = map.get(seqId);
      } else {
        st = new TreeSet<>();
        map.put(seqId, st);
      }
      final int start = Integer.parseInt(sStart);
      final int end = Integer.parseInt(sEnd);
      st.add(new SimpleCnvRegion(start, end));
    }
    final Map<String, Region> mapr = new HashMap<>();
    for (final Map.Entry<String, SortedSet<AbstractCnvRegion>> entry : map.entrySet()) {
      final SortedSet<AbstractCnvRegion> set = entry.getValue();
      final int size = set.size();
      assert size != 0;
      final int minStart = set.first().getStart();
      int maxEnd = Integer.MIN_VALUE;
      for (final AbstractCnvRegion region : set) {
        maxEnd = Math.max(maxEnd, region.getEnd());
      }
      final Region res;
      if (size == 1) {
        res = new SimpleCnvRegion(minStart, maxEnd);
      } else {
        res = new ComplexCnvRegion(minStart, maxEnd, set);
      }
      mapr.put(entry.getKey(), res);
    }
    return mapr;
  }

  /**
   * Method for producing N-Region map from Sequences.
   * @param reader the Sequences reader
   * @param blockSize minimum block size for N-Regions
   * @return the region object (never null)
   * @throws IOException when reader fails
   * @throws IllegalStateException when not a DNA sequence
   */
  public static Map<String, Region> regionsFromSDF(final SequencesReader reader, int blockSize) throws IOException {
    if (reader.type() == SequenceType.DNA) {
      final Map<String, Region> map = new HashMap<>();
      final byte[] buf = new byte[(int) reader.maxLength()];
      for (long i = 0; i < reader.numberSequences(); ++i) {
        final int length = reader.read(i, buf);
        map.put(reader.name(i), detectNs(buf, 0, length, blockSize));
      }
      return map;
    } else {
      throw new UnsupportedOperationException("Only support DNA sequences");
    }
  }

  static Region detectNs(byte[] seq, int start, int end, int blockSize) {
    final SortedSet<AbstractCnvRegion> set = new TreeSet<>();
    int nStart = 0;
    boolean nRegion = false;
    for (int i = start; i < end; ++i) {
      final boolean n = seq[i] <= MAX_UNKNOWN;
      if (!nRegion && n) {
        nStart = i;
        nRegion = true;
      } else if (nRegion && !n) {
        if ((i - nStart) >= blockSize) {
          set.add(new SimpleCnvRegion(nStart / blockSize + 1, (i - (nStart % blockSize)) / blockSize + 1));
        }
        nRegion = false;
      }
    }
    if (nRegion) {
      if ((end - nStart) >= blockSize) {
        set.add(new SimpleCnvRegion(nStart / blockSize + 1, (end - (nStart % blockSize)) / blockSize + 1));
      }
    }
    final Region res;
    if (set.size() == 0) {
      res = EmptyRegion.EMPTY_REGION;
    } else if (set.size() == 1) {
      res = set.first();
    } else {
      res = new ComplexCnvRegion(set.first().getStart(), set.last().getEnd(), set);
    }
    return res;
  }

  /**
   * Calculate the germ-line deletions under the factor of the mean
   * @param bucketCounts the bucket count numbers for the germ-line
   * @param factor the factor to divide the mean by
   * @param ignore region to ignore
   * @param penaltyOn switch off adding extra penalty
   * @return Region defining the germ-line deletions
   */
  public static Region findGermlineDeletesUnderMean(int[] bucketCounts, double factor, Region ignore, boolean penaltyOn) {
    final double meandiv = ArrayUtils.sum(bucketCounts) / (double) bucketCounts.length / factor;
    final SortedSet<AbstractCnvRegion> set = new TreeSet<>();
    final LeastSquaresModel lsm = new LeastSquaresModel("", 2, 1, 3.0, penaltyOn, new LsmOutput() {
      @Override
      public void close() { }
      @Override
      public void out(String id, int start, int end, double sum, double sum2) {
        final double mean = sum / (double) (end - start);
        if (mean < meandiv) {
          set.add(new SimpleCnvRegion(start, end));
        }
      }
    });
    return findGermlineDeletes(bucketCounts, ignore, set, lsm);
  }

  /**
   * Calculate the germ-line deletions over the factor of the mean
   * @param bucketCounts the bucket count numbers for the germ-line
   * @param factor the factor to multiply the mean by
   * @param ignore region to ignore
   * @param penaltyOn switch penalty on
   * @return Region defining the germ-line deletions
   */
  public static Region findGermlineDeletesOverMean(int[] bucketCounts, double factor, Region ignore, boolean penaltyOn) {
    final double meanmul = factor * (ArrayUtils.sum(bucketCounts) / (double) bucketCounts.length);
    final SortedSet<AbstractCnvRegion> set = new TreeSet<>();
    final LeastSquaresModel lsm = new LeastSquaresModel("", 2, 1, 3.0, penaltyOn, new LsmOutput() {
      @Override
      public void close() { }
      @Override
      public void out(String id, int start, int end, double sum, double sum2) {
        final double mean = sum / (double) (end - start);
        if (mean > meanmul) {
          set.add(new SimpleCnvRegion(start, end));
        }
      }
    });
    return findGermlineDeletes(bucketCounts, ignore, set, lsm);
  }

  private static Region findGermlineDeletes(int[] bucketCounts, Region ignore, final SortedSet<AbstractCnvRegion> set, final LeastSquaresModel lsm) {
    assert lsm != null;
    int blocks = 0;
    int seqStart = -1;
    State state = State.OUT;
    for (final int bucket : bucketCounts) {
      final boolean okRegion = !ignore.isInRegion(blocks + 1);
      if (state == State.IN) {
        if (!okRegion) {
          try {
            lsm.scan(seqStart, blocks + 1);
          } catch (IOException e) {
            //Unpossible
          }
          seqStart = -1;
          state = State.OUT;
        }
      } else {
        if (okRegion) {
          seqStart = blocks + 1;
          state = State.IN;
        }
      }
      lsm.add((double) bucket);
      ++blocks;
    }
    if (state == State.IN) {
      try {
        lsm.scan(seqStart, blocks + 1);
      } catch (IOException e) {
        //Unpossible
      }
    }
    final Region res;
    if (set.size() == 0) {
      res = EmptyRegion.EMPTY_REGION;
    } else if (set.size() == 1) {
      res = set.first();
    } else {
      res = new ComplexCnvRegion(set.first().getStart(), set.last().getEnd(), set);
    }
    return res;
  }
}
