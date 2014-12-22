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
import com.rtg.util.integrity.Exam;
import com.rtg.util.integrity.Integrity;

/**
 * A wrapper for a sequences reader which reverse complements bases and reverse
 * quality information.
 */
public final class ReverseComplementingReader implements SequencesReader, Integrity {

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

  @Override
  public long maxLength() {
    return mUnderlyingReader.maxLength();
  }

  @Override
  public long minLength() {
    return mUnderlyingReader.minLength();
  }

  @Override
  public boolean hasQualityData() {
    return mUnderlyingReader.hasQualityData();
  }

  @Override
  public boolean hasNames() {
    return mUnderlyingReader.hasNames();
  }

  @Override
  public int readCurrent(final byte[] dataOut) throws IllegalArgumentException, IOException {
    final int r = mUnderlyingReader.readCurrent(dataOut);
    DNA.reverseComplementInPlace(dataOut, 0, r);
    return r;
  }

  @Override
  public int readCurrent(byte[] dataOut, int start, int length) throws IllegalArgumentException, IOException {
    final int r = mUnderlyingReader.readCurrent(dataOut, start, length);
    DNA.reverseComplementInPlace(dataOut, start, start + r);
    return r;
  }

  @Override
  public void seek(final long sequenceId) throws IOException {
    mUnderlyingReader.seek(sequenceId);
  }

  @Override
  public boolean nextSequence() throws IOException {
    return mUnderlyingReader.nextSequence();
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
  public int readCurrentQuality(final byte[] dest) throws IllegalArgumentException, IllegalStateException, IOException {
    final int r = mUnderlyingReader.readCurrentQuality(dest);
    reverse(dest, 0, r);
    return r;
  }

  @Override
  public int readCurrentQuality(byte[] dest, int start, int length) throws IllegalArgumentException, IllegalStateException, IOException {
    final int r = mUnderlyingReader.readCurrentQuality(dest, start, length);
    reverse(dest, start, length);
    return r;
  }

  @Override
  public int currentLength() {
    return mUnderlyingReader.currentLength();
  }

  @Override
  public String currentName() throws IOException {
    return mUnderlyingReader.currentName();
  }

  @Override
  public String currentFullName() throws IOException {
    return currentName();
  }

  @Override
  public long currentSequenceId() {
    return mUnderlyingReader.currentSequenceId();
  }

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
  public int length(final long sequenceIndex) {
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
  public long numberSequences() {
    return mUnderlyingReader.numberSequences();
  }

  @Override
  public long totalLength() {
    return mUnderlyingReader.totalLength();
  }

  @Override
  public SequenceType type() {
    return mUnderlyingReader.type();
  }

  @Override
  public boolean integrity() {
    return Exam.integrity(mUnderlyingReader);
  }

  @Override
  public boolean globalIntegrity() {
    return integrity();
  }

  @Override
  public File path() {
    return mUnderlyingReader.path();
  }

  @Override
  public void close() throws IOException {
    mUnderlyingReader.close();
  }

  private void swap(final long[] x, final int a, final int b) {
    final long t = x[a];
    x[a] = x[b];
    x[b] = t;
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
  public long dataChecksum() {
    return mUnderlyingReader.dataChecksum();
  }
  @Override
  public long qualityChecksum() {
    return mUnderlyingReader.qualityChecksum();
  }
  @Override
  public long nameChecksum() {
    return mUnderlyingReader.nameChecksum();
  }

  @Override
  public PrereadNamesInterface names() throws IOException {
    return mUnderlyingReader.names();
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
  public long[] histogram() {
    return mUnderlyingReader.histogram();
  }

  @Override
  public long[] posHistogram() {
    return mUnderlyingReader.posHistogram();
  }

  @Override
  public double globalQualityAverage() {
    return mUnderlyingReader.globalQualityAverage();
  }

  @Override
  public boolean hasHistogram() {
    return mUnderlyingReader.hasHistogram();
  }

  @Override
  public long longestNBlock() {
    return mUnderlyingReader.longestNBlock();
  }

  @Override
  public long nBlockCount() {
    return mUnderlyingReader.nBlockCount();
  }

  @Override
  public PrereadArm getArm() {
    return mUnderlyingReader.getArm();
  }

  @Override
  public PrereadType getPrereadType() {
    return mUnderlyingReader.getPrereadType();
  }

  @Override
  public SdfId getSdfId() {
    return mUnderlyingReader.getSdfId();
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
  public long sdfVersion() {
    return mUnderlyingReader.sdfVersion();
  }

  @Override
  public SequencesReader copy() {
    return new ReverseComplementingReader(mUnderlyingReader.copy());
  }

  @Override
  public boolean compressed() {
    return mUnderlyingReader.compressed();
  }

  @Override
  public void reset() {
    mUnderlyingReader.reset();
  }

  @Override
  public String currentNameSuffix() throws IllegalStateException, IOException {
    return mUnderlyingReader.currentNameSuffix();
  }

  @Override
  public String nameSuffix(long sequenceIndex) throws IOException {
    return mUnderlyingReader.nameSuffix(sequenceIndex);
  }

  @Override
  public long suffixChecksum() {
    return mUnderlyingReader.suffixChecksum();
  }

  @Override
  public String getReadMe() throws IOException {
    return mUnderlyingReader.getReadMe();
  }

}

