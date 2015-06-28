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
package com.rtg.simulation.genome;

import java.io.File;
import java.io.IOException;

import com.rtg.mode.DNA;
import com.rtg.mode.Residue;
import com.rtg.mode.SequenceType;
import com.rtg.reader.PrereadType;
import com.rtg.reader.SequenceDataSource;
import com.rtg.reader.SequencesWriter;
import com.rtg.util.Constants;
import com.rtg.util.PortableRandom;

/**
 * Generates a randomized DNA sequence
 *
 *
 */
public class SequenceGenerator {

  private final SequencesWriter mWriter;

  private static class RandomDataSource implements SequenceDataSource {

    private int mSequenceNumber = 0;
    private final int[] mLengths;
    private final PortableRandom mSource;
    private final RandomDistribution mDistribution;
    private long mMaxLength = Long.MIN_VALUE;
    private long mMinLength = Long.MAX_VALUE;
    private final Residue[] mResidues = DNA.values();
    private final String mPrefix;

    private byte[] mSequenceData;

    public RandomDataSource(final int[] lengths, final PortableRandom source, final RandomDistribution distribution, String namePrefix) {
      mLengths = lengths;
      mSource = source;
      mDistribution = distribution;
      mPrefix = namePrefix;
    }

    @Override
    public void close() {
    }

    @Override
    public String name() {
      return mPrefix + mSequenceNumber;
    }

    @Override
    public SequenceType type() {
      return SequenceType.DNA;
    }

    @Override
    public long getWarningCount() {
      return 0;
    }

    @Override
    public boolean hasQualityData() {
      return false;
    }

    @Override
    public boolean nextSequence() {
      if (mSequenceNumber >= mLengths.length) {
        return false;
      }
      mSequenceData = new byte[mLengths[mSequenceNumber]];
      for (int i = 0; i < mSequenceData.length; i++) {
        mSequenceData[i] = (byte) getRandomResidue().ordinal();
      }

      mSequenceNumber++;
      mMinLength = Math.min(mMinLength, currentLength());
      mMaxLength = Math.max(mMaxLength, currentLength());
      return true;
    }

    @Override
    public byte[] sequenceData() {
      return mSequenceData;
    }


    private Residue getRandomResidue() {
      if (mDistribution == null) {
        return mResidues[1 + mSource.nextInt(mResidues.length - 1)];
      }
      return mResidues[1 + mDistribution.nextValue()];
    }

    @Override
    public byte[] qualityData() {
      return null;
    }

    @Override
    public void setDusting(final boolean val) { }

    @Override
    public int currentLength() {
      return mSequenceData.length;
    }

    @Override
    public long getDusted() {
      return 0;
    }

    @Override
    public long getMaxLength() {
      return mMaxLength;
    }

    @Override
    public long getMinLength() {
      return mMinLength;
    }
  }

  public long getSizeLimit() {
    return mWriter.getSizeLimit();
  }

  /**
   * Create a sequence generator with specified seed and nucleotide frequency
   * corresponding to the supplied distribution
   *  @param generator random number generator
   * @param distribution relative frequency of a/t/g/c
   * @param lengths array of lengths of the sequences to write
   * @param outDirectory output directory
   * @param namePrefix prefix to use for generated sequence names
   */
  public SequenceGenerator(final PortableRandom generator, final RandomDistribution distribution, final int[] lengths, final File outDirectory, String namePrefix) {
    final RandomDataSource dataSource = new RandomDataSource(lengths, generator, distribution, namePrefix);
    mWriter = new SequencesWriter(dataSource, outDirectory, Constants.MAX_FILE_SIZE, PrereadType.UNKNOWN, true);
  }

  /**
   * @param comment the comment for the generated SDF file(s)
   */
  public void setComment(final String comment) {
    mWriter.setComment(comment);
  }

  /**
   * Write a sequences
   *
    * @throws IOException when IO errors occur
   */
  public void createSequences() throws IOException {
    mWriter.processSequences();
  }

}
