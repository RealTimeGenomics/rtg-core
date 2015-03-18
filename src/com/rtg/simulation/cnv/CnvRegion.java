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
package com.rtg.simulation.cnv;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Locale;

import com.rtg.reader.SequencesReader;
import com.rtg.util.PortableRandom;
import com.rtg.util.StringUtils;
import com.rtg.util.Utils;

/**
 */
public class CnvRegion {

  private static final String TB = "\t";

  boolean mIsCnved;
  final int mSequenceId;
  final int mStart;
  int mLength;
  int mNumCopies;
  boolean mOutputDelete;
  boolean mTwinDelete;
  final ArrayList<CnvRegion> mCopies;
  final ArrayList<Boolean> mCopiesHetero;
  CnvPriorParams mPriors;

  CnvRegion(final int sequenceId, final int start, int length, CnvPriorParams priors) {
    mIsCnved = false;
    mSequenceId = sequenceId;
    mStart = start;
    mLength = length;
    mPriors = priors;
    mOutputDelete = false;
    mTwinDelete = false;

    mCopies = new ArrayList<>();
    mCopiesHetero = new ArrayList<>();
    mNumCopies = 0;
    //System.out.println("sequence=" + mSequenceId + " start=" + mStart + " length=" + mLength);
  }

  /**
   * add a region to merge
   * @param testedRegion region to be added
   */
  void addRegion(CnvRegion testedRegion) {
    mLength += testedRegion.mLength;
  }

  boolean isUnpickable(int[] sequenceLengths, boolean dontIgnoreCnved, int minLength) {
    // dont pick if already picked
    // or has length zero = is end region
    // of is a region over the whole sequence
    if ((dontIgnoreCnved && mIsCnved) || mLength == 0
        || mLength == sequenceLengths[mSequenceId]) {
      return true;
    }
    // test length again for cuts at end of sequence
    // test if allowed per priors
    if (mPriors != null) {
      final int index = CnvSimulator.getMagnitudeIndex(mLength, mPriors.powerLengthDistribution().length - 1);
      if (mPriors.powerLengthDistribution()[index] == 0.0) {
        return true;
      }
    } else if (mLength < minLength) {
      // test if allowed per min length
      return true;
    }
    return false;
  }

  void initializeAsCnved(PortableRandom random) {
    mIsCnved = true;
    final double rand = random.nextDouble();
    if (rand < mPriors.probDeletedOnBothStrands()) {
      mOutputDelete = true;
      mTwinDelete = true;
    } else if (rand < mPriors.probDeletedOnBothStrands() + mPriors.probDeletedOnOneStrand()) {
      if (random.nextDouble() < 0.5) {
        mOutputDelete = true;
        mTwinDelete = false;
      } else {
        mTwinDelete = true;
        mOutputDelete = false;
      }
    } else {
      mNumCopies = generateNumCopies(random, mPriors, mLength);
      mOutputDelete = false;
      mTwinDelete = false;
    }
  }

  static int generateNumCopies(PortableRandom random, CnvPriorParams priors, int length) {
    final int magnitudeIndex = CnvSimulator.getMagnitudeIndex(length, 7);
    final double[] distribution = priors.copyNumberThresholds()[magnitudeIndex];
    int index;
    index = CnvSimulator.chooseFromAccumDistribution(random.nextDouble(), distribution);
    if (index == 0) {
      return 1;
    } else {
      return random.nextInt(priors.copyRangeEnds()[index] - priors.copyRangeEnds()[index - 1])
      + priors.copyRangeEnds()[index - 1] + 1;
    }
  }

  /**
   * @param k index of the copy to get
   * @param input SDF file to read from
   * @return bytes to be copied
   * @throws IOException if error in IO
   * @throws IllegalStateException if illegal state
   * @throws IllegalArgumentException if illegal argument
   */
  byte[] getCopy(int k, SequencesReader input) throws IllegalArgumentException, IllegalStateException, IOException {
    final CnvRegion copy = mCopies.get(k);
    final byte[] genomeSeq = new byte[input.length(copy.mSequenceId)];
    input.read(copy.mSequenceId, genomeSeq);
    return CnvSimulator.getDnaArray(genomeSeq, copy.mStart, copy.mLength);
  }

  /**
   * @param k index of the copy to get
   * @return if the copy is done on both strands
   */
  boolean getBothStrands(int k) {
    return mCopiesHetero.get(k);
  }

  static final double ERROR = 0.0;
  static final int DB_CN = 0;
  private static final String LS = StringUtils.LS;

  CnvRegion(final int sequenceId, final int start, int length) {
    mSequenceId = sequenceId;
    mStart = start;
    mLength = length;
    mOutputDelete = false;
    mTwinDelete = false;
    mNumCopies = 0;
    mCopies = new ArrayList<>();
    mCopiesHetero = new ArrayList<>();

    //System.out.println("sequence=" + mSequenceId + " start=" + mStart + " length=" + mLength);
  }

  CnvRegion(final int sequenceId, final int start, int length, int copyNumber) {
    mIsCnved = true;
    mSequenceId = sequenceId;
    mStart = start;
    mLength = length;
    if (copyNumber < 2) {
      mNumCopies = 0;
      mTwinDelete = true;
      mOutputDelete = copyNumber < 1;
    } else {
      mNumCopies = copyNumber - 2;
    }
    mCopies = new ArrayList<>();
    mCopiesHetero = new ArrayList<>();

    //System.out.println("sequence=" + mSequenceId + " start=" + mStart + " length=" + mLength);
  }

  /**
   * @param currentRegion current region
   * @param random the random number generator to use to order the copies
   * @param isHetero is heterozygous
   */
  void addCopy(CnvRegion currentRegion, PortableRandom random, boolean isHetero) {
    final int index = random.nextInt(mCopies.size() + 1);
    mCopies.add(index, currentRegion);
    mCopiesHetero.add(index, isHetero);

  }

  /**
   * @param currentRegion the current region
   */
  void addCopy(CnvRegion currentRegion) {
    mCopies.add(currentRegion);
    mCopiesHetero.add(false);
  }

  /**
   * @return the copy number for this region
   */
  protected int getCN() {
    if (mLength == 0) {
      return 0;
    }
    int num = 0;
    if (!mOutputDelete) {
      num++;
    }
    if (!mTwinDelete) {
      num++;
    }
    num += mNumCopies;
    return num;
  }

  /**
   * @param seqId sequence id
   * @param seqName sequence name
   * @return a byte array containing a string representation of this region
   */
  public byte[] toBytes(int seqId, String seqName) {
    assert seqId == mSequenceId;
    final StringBuilder line = new StringBuilder();
    line.append(seqName).append(TB);
    line.append(mStart).append(TB);
    line.append(mStart + mLength).append(TB);
    line.append("cnv").append(TB);
    line.append(getCN()).append(TB);
    line.append(DB_CN).append(TB);
    line.append(Utils.realFormat(ERROR)).append(TB).append(LS);
    return line.toString().getBytes();
  }

  /**
   * @return a byte array containing a string representation of this region
   */
  @Override
  public String toString() {
    final StringBuilder line = new StringBuilder();
    line.append("seq: ").append(mSequenceId);
    line.append(" cnved: ").append(mIsCnved);
    line.append(" start: ").append(mStart);
    line.append(" length: ").append(mLength);
    line.append(" del-1: ").append(mOutputDelete);
    line.append(" del-2: ").append(mTwinDelete);
    line.append(" copies: ").append(mNumCopies);
    line.append(" num copies added: ").append(mCopies.size());
    line.append(" num flags added: ").append(mCopiesHetero.size()).append(LS);
    line.append(" HeteroList:");
    for (Boolean aMCopiesHetero : mCopiesHetero) {
      line.append(" ").append(aMCopiesHetero);
    }
    return line.toString().toLowerCase(Locale.getDefault());
  }


}
