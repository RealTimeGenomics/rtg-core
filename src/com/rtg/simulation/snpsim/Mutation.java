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
package com.rtg.simulation.snpsim;

import java.util.Arrays;

import com.rtg.reference.Ploidy;
import com.rtg.util.EnumHelper;
import com.rtg.util.PortableRandom;
import com.rtg.util.PseudoEnum;

/**
 * Defines a genome mutation.
 *
 */
public class Mutation {

  final int mPos;
  final MutationType mType;
  final boolean mHeterozygous;
  final Ploidy mPloidy;
  final DifferentMode mDiffMode; // not used after construction
  final GenDiffMode mGenDiffMode;
  final int mLength;
  final int mLengthTwin;
  final boolean mFirstDelete; // only for insertion-deletions
  final byte[] mBases;

  /**
   * Constructor for a genome mutation
   *
   * @param r random number generator
   * @param pos position of the mutation
   * @param gen proxy to parameter priors
   * @param sequenceLength length of current sequence
   * @param refBase the base at the position to be generate, is used only when generate SNPs
   * @param ploidy the ploidy of the mutation location
   */
  public Mutation(PortableRandom r, int pos, GenomeMutatorPriors gen, int sequenceLength, byte refBase, Ploidy ploidy) {
    this(r, pos, gen, sequenceLength, refBase, gen.chooseType(r), ploidy);
  }

  /**
   * @param r the base
   * @param pos the position of the mutation
   * @param gen the priors to use to generate the mutation
   * @param sequenceLength the length of the current sequence
   * @param refBase the bases at the reference position
   * @param chooseType the type of the mutation
   * @param ploidy the ploidy of the mutation location
   */
  public Mutation(PortableRandom r, int pos, GenomeMutatorPriors gen, int sequenceLength, byte refBase, MutationType chooseType, Ploidy ploidy) {
    mType = chooseType;
    mPos = pos;
    mPloidy = ploidy;
    boolean firstDelete = true;
    if (mType == MutationType.SNP && gen.hasPriors()) {
      final byte[] bases = gen.chooseAltSnp(r, refBase);
      if (bases[0] == bases[1]) {
        mHeterozygous = false;
        mDiffMode = DifferentMode.HOMOZYGOUS;
        mGenDiffMode = GenDiffMode.BOTH_SAME;
        mBases = bases;
      } else {
        mHeterozygous = true;

        if (refBase == bases[0] || refBase == bases[1]) {
          mDiffMode = DifferentMode.ONE_ONLY;
          if (r.nextBoolean()) {
            mGenDiffMode = GenDiffMode.FIRST_ONLY;
            if (refBase == bases[0]) {
              mBases = new byte[] {bases[1], bases[0]};
            } else {
              mBases = new byte[] {bases[0], bases[1]};
            }
          } else {
            mGenDiffMode = GenDiffMode.TWIN_ONLY;
            if (refBase == bases[0]) {
              mBases = new byte[] {bases[1], bases[0]};
            } else {
              mBases = new byte[] {bases[0], bases[1]};
            }
          }
        } else {
          if (r.nextBoolean()) {
            mBases = new byte[] {bases[0], bases[1]};
          } else {
            mBases = new byte[] {bases[1], bases[0]};
          }
          mDiffMode = DifferentMode.DIFFERENT;
          mGenDiffMode = GenDiffMode.DIFFERENT;
        }
      }
    } else {
      mBases = null;
      DifferentMode selectedDiffMode = gen.chooseDifferentMode(r, mType);

      // Heterozygous deletes at the end of the sequence are tricky.
      // Heterozygous both different makes no sense in that case
      while (sequenceLength - pos == 1 && mType == MutationType.DELETE && selectedDiffMode != DifferentMode.ONE_ONLY && selectedDiffMode != DifferentMode.HOMOZYGOUS) {
        selectedDiffMode = gen.chooseDifferentMode(r, mType);
      }
      mDiffMode = selectedDiffMode;
      mHeterozygous = mDiffMode != DifferentMode.HOMOZYGOUS;
      if (!mHeterozygous) {
        mGenDiffMode = GenDiffMode.BOTH_SAME;
      } else if (mDiffMode == DifferentMode.ONE_ONLY) {
        if (r.nextBoolean()) {
          mGenDiffMode = GenDiffMode.FIRST_ONLY;
        } else {
          mGenDiffMode = GenDiffMode.TWIN_ONLY;
        }
      } else {
        mGenDiffMode = GenDiffMode.DIFFERENT;
      }
    }
    // need the lengths immediately to be able to check for overlaps
    int length1 = 0;
    int length2 = 0;
    switch (mType.simple()) {
      case SNP:
        length1 = 1;
        length2 = 1;
        break;
      case MNP:
      case INSERT:
      case DELETE:
      case INSDEL:
        firstDelete = r.nextBoolean();
        final int len = gen.chooseLength(r, mType, mHeterozygous);
        if (!mHeterozygous) {
          length1 = len;
          length2 = len;
        } else {
          switch (mGenDiffMode) {
            case FIRST_ONLY:
              length1 = len;
              break;
            case TWIN_ONLY:
              length2 = len;
              break;
            case DIFFERENT:
              length1 = len;
              length2 = gen.chooseLength(r, mType, mHeterozygous);
              break;
            default:
              break;
          }
        }
        break;
      default:
        throw new IllegalStateException();
    }
    if (mType != MutationType.INSERT) {
      if (length1 > sequenceLength - pos) {
        length1 = sequenceLength - pos;
      }
      if (length2 > sequenceLength - pos) {
        length2 = sequenceLength - pos;
      }
      if (mType == MutationType.DELETE && mHeterozygous && length1 == length2) {
        if (r.nextBoolean()) {
          length2--;
        } else {
          length1--;
        }
      }
    }
    // System.out.println(mType + " " + length1 + " " + length2);
    mLength = length1;
    mLengthTwin = length2;
    mFirstDelete = firstDelete;
  }

  /**
   * Constructor with all fields as input parameter
   *  @param pos position of mutation
   * @param type type of mutation
   * @param heterozygous if heterozygous
   * @param diffMode mode if both variations of are different
   * @param genDiffMode modes for generation
   * @param length length of first variation
   * @param length2 length of second variation
   * @param bases the bases created using prior based SNP generation
   * @param ploidy the ploidy of the mutation location
   */
  public Mutation(int pos, MutationType type, boolean heterozygous, DifferentMode diffMode, GenDiffMode genDiffMode, int length, int length2, byte[] bases, Ploidy ploidy) {
    mPos = pos;
    mType = type;
    mPloidy = ploidy;
    mHeterozygous = heterozygous;
    mDiffMode = diffMode;
    mGenDiffMode = genDiffMode;
    mLength = length;
    mLengthTwin = length2;
    mFirstDelete = true;
    if (bases != null) {
      mBases = Arrays.copyOf(bases, 2);
    } else {
      mBases = null;
    }
  }

  /** Enumeration of mutation types. */
  public enum SimpleMutationType {
    /** SNP mutation */
    SNP,
    /** multiple SNP mutation */
    MNP,
    /** insertion */
    INSERT,
    /** deletion */
    DELETE,
    /** heterozygous insert plus deletion */
    INSDEL,
    /** mutation read from SNP file */
    PREPARED,
  }

  /**
   * Class for mutation types that are generated
   */
  public static final class MutationType implements PseudoEnum {
    /** SNP mutation */
    public static final MutationType SNP = new MutationType(SimpleMutationType.SNP);
    /** multiple SNP mutation */
    public static final MutationType MNP = new MutationType(SimpleMutationType.MNP);
    /** insertion */
    public static final MutationType INSERT = new MutationType(SimpleMutationType.INSERT);
    /** deletion */
    public static final MutationType DELETE = new MutationType(SimpleMutationType.DELETE);
    /** heterozygous insert plus deletion */
    public static final MutationType INSDEL = new MutationType(SimpleMutationType.INSDEL);
    /** mutation read from SNP file */
    public static final MutationType PREPARED = new MutationType(SimpleMutationType.PREPARED);

    private final SimpleMutationType mBase;

    private static final EnumHelper<MutationType> HELPER = new EnumHelper<>(MutationType.class, new MutationType[] {SNP, MNP, INSERT, DELETE, INSDEL, PREPARED});

    private MutationType(SimpleMutationType base) {
      mBase = base;
    }

    /**
     * @return list of the enum names
     */
    public static String[] names() {
      return HELPER.names();
    }

    @Override
    public String name() {
      return mBase.toString();
    }

    @Override
    public String toString() {
      return mBase.toString();
    }

    @Override
    public int ordinal() {
      return mBase.ordinal();
    }

    /**
     * Return switchable version
     *
     * @return the base version
     */
    public SimpleMutationType simple() {
      return mBase;
    }

    /**
     * {@link EnumHelper#valueOf(java.lang.String)}
     *
     * @param str as in link
     * @return as in link
     */
    public static MutationType valueOf(String str) {
      return HELPER.valueOf(str);
    }

    /**
     * {@link EnumHelper#values()}
     *
     * @return as in link
     */
    public static MutationType[] values() {
      return HELPER.values();
    }
  }

  /** Enumeration of modes. */
  public enum DifferentMode {
    /** homozygous means both are the same at this mutation point */
    HOMOZYGOUS,
    /** heterozygous and only one genome is different */
    ONE_ONLY,
    /** heterozygous and both genomes are different */
    DIFFERENT,
  }

  String makeString() {
    final StringBuilder str = new StringBuilder();
    str.append("   pos:").append(mPos).append(" ");
    if (mHeterozygous) {
      str.append("e ");
    } else {
      str.append("o ");
    }
    str.append(mType.toString()).append(" ");
    str.append("l1: ").append(mLength).append(" l2: ").append(mLengthTwin).append(" ");
    str.append(mGenDiffMode.toString());

    return str.toString();
  }


  boolean isFirstDelete() {
    return mFirstDelete;
  }

}
