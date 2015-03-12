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
package com.rtg.reader;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;

import com.rtg.mode.DNA;
import com.rtg.mode.SequenceType;

/**
 * A wrapper for a sequences reader which reverse complements bases and reverse
 * quality information.
 */
public final class ReverseComplementingReader extends AbstractSequencesReader {

  private final SequencesReader mUnderlyingReader;

  /**
   * Construct a new reverse complementing reader wrapping another reader.
   *
   * @param reader a reader
   * @exception NullPointerException if <code>reader</code> is null
   * @exception IllegalArgumentException if the underlying type is not nucleotides
   */
  public ReverseComplementingReader(final SequencesReader reader) {
    if (reader == null) {
      throw new NullPointerException();
    }
    if (reader.type() != SequenceType.DNA) {
      throw new IllegalArgumentException();
    }
    mUnderlyingReader = reader;
  }


  static void reverse(final byte[] x, final int start, final int length) {
    int l = start;
    int r = start + length - 1;
    while (l < r) {
      final byte b = x[l];
      x[l] = x[r];
      x[r] = b;
      l++;
      r--;
    }
  }

  @Override
  public IndexFile index() {
    return mUnderlyingReader.index();
  }


  // Direct accessor methods
  @Override
  public int readQuality(final long sequenceIndex, final byte[] dest) throws IllegalArgumentException, IOException {
    final int r = mUnderlyingReader.readQuality(sequenceIndex, dest);
    reverse(dest, 0, r);
    return r;
  }

  @Override
  public int readQuality(long sequenceIndex, byte[] dest, int start, int length) throws IllegalArgumentException, IOException {
    final int r = mUnderlyingReader.readQuality(sequenceIndex, dest, start, length);
    reverse(dest, 0, r);
    return r;
  }

  @Override
  public int read(final long sequenceIndex, final byte[] dataOut) throws IllegalArgumentException, IOException {
    final int r = mUnderlyingReader.read(sequenceIndex, dataOut);
    DNA.reverseComplementInPlace(dataOut, 0, r);
    return r;
  }

  @Override
  public int read(final long sequenceIndex, final byte[] dataOut, int start, int length) throws IllegalArgumentException, IOException {
    final int r = mUnderlyingReader.read(sequenceIndex, dataOut, start, length);
    DNA.reverseComplementInPlace(dataOut, start, start + r);
    return r;
  }

  @Override
  public int length(final long sequenceIndex) throws IOException {
    return mUnderlyingReader.length(sequenceIndex);
  }

  @Override
  public String name(final long sequenceIndex) throws IOException {
    return mUnderlyingReader.name(sequenceIndex);
  }

  @Override
  public String fullName(final long sequenceIndex) throws IOException {
    return name(sequenceIndex);
  }

  @Override
  public File path() {
    return mUnderlyingReader.path();
  }

  @Override
  public void close() throws IOException {
    mUnderlyingReader.close();
  }

  @Override
  public PrereadNamesInterface names() throws IOException {
    return mUnderlyingReader.names();
  }

  private void swap(final long[] x, final int a, final int b) {
    final long t = x[a];
    x[a] = x[b];
    x[b] = t;
  }

  @Override
  public long numberSequences() {
    return mUnderlyingReader.numberSequences();
  }

  @Override
  public long[] residueCounts() {
    final long[] q = mUnderlyingReader.residueCounts();
    final long[] c = Arrays.copyOf(q, q.length);
    swap(c, DNA.A.ordinal(), DNA.T.ordinal());
    swap(c, DNA.C.ordinal(), DNA.G.ordinal());
    return c;
  }

  @Override
  public long lengthBetween(final long start, final long end) throws IOException {
    return mUnderlyingReader.lengthBetween(start, end);
  }

  /**
   * Puts lengths of sequences in an array and returns it.
   * Lightly tested
   *
   * @param start starting sequence
   * @param end ending sequence (excl)
   * @return array of lengths
   * @throws IOException if an I/O error occurs
   */
  @Override
  public int[] sequenceLengths(final long start, final long end) throws IOException {
    return mUnderlyingReader.sequenceLengths(start, end);
  }

  @Override
  public double[] positionQualityAverage() {
    final double[] av = mUnderlyingReader.positionQualityAverage();
    int l = 0;
    int r = av.length - 1;
    while (l < r) {
      final double b = av[l];
      av[l] = av[r];
      av[r] = b;
      l++;
      r--;
    }
    return av;
  }

  @Override
  public SequencesReader copy() {
    return new ReverseComplementingReader(mUnderlyingReader.copy());
  }

  @Override
  public String nameSuffix(long sequenceIndex) throws IOException {
    return mUnderlyingReader.nameSuffix(sequenceIndex);
  }

  @Override
  public String getReadMe() throws IOException {
    return mUnderlyingReader.getReadMe();
  }

}

