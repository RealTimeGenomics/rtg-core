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
package com.rtg.variant;

import java.util.Arrays;

import com.rtg.launcher.globals.CoreGlobalFlags;
import com.rtg.launcher.globals.GlobalFlags;
import com.rtg.mode.DNA;
import com.rtg.mode.DnaUtils;
import com.rtg.reader.Arm;
import com.rtg.reader.CgUtils;
import com.rtg.reader.FastaUtils;
import com.rtg.sam.MateInfo;
import com.rtg.sam.ReaderRecord;
import com.rtg.sam.SamUtils;
import com.rtg.util.CompareHelper;
import com.rtg.util.MathUtils;
import com.rtg.util.intervals.SequenceIdLocusSimple;

import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;

/**
 * Alignment information needed by variant implementations.
 */
public final class VariantAlignmentRecord extends SequenceIdLocusSimple implements ReaderRecord<VariantAlignmentRecord>, MateInfo, MapInfo {
  // Factor 1.1 to cover arithmetic error in sum used for error() in EvidenceComplex
  private static final boolean MIN_QUALITY_AS_TWO = GlobalFlags.getBooleanValue(CoreGlobalFlags.MIN_BASE_QUALITY_AS_TWO);

  private static final int FLAG_MATED = 1;
  private static final int FLAG_PAIRED = 2;
  private static final int FLAG_FIRST = 4;
  private static final int FLAG_NEGATIVE = 8;
  private static final int FLAG_UNMAPPED = 16;

  /**
   * Record used to denote an overflow condition, not a true record.
   * @param start 0-based start position of overflow
   * @param length of overflow
   * @return overflow record
   */
  public static VariantAlignmentRecord overflow(final int start, final int length) {
    return new VariantAlignmentRecord(start, start + length);
  }

  private final byte[] mBases;
  private final byte[] mRecalibratedQuality;
  private final String mCigar;
  private final byte mMappingQuality;
  private final byte mFlag; // Not the same semantics as SAM flag.
  private final int mAmbiguity;
  private final int mAlignmentScore;
  private final int mMateSequenceId;
  private final int mFragmentLength;

  private final SAMReadGroupRecord mReadGroup;
  private final String mSuperCigar;
  private final byte[] mOverlapBases;
  private final byte[] mOverlapQuality;
  private final String mOverlapInstructions;
  private final String mCgReadDelta;

  private final int mGenome;

  private VariantAlignmentRecord(final int start, final int end) {
    super(0, start, end);
    mBases = null;
    mRecalibratedQuality = null;
    mCigar = null;
    mMappingQuality = 0;
    mFlag = -1;
    mAmbiguity = 0;
    mAlignmentScore = 0;
    mMateSequenceId = 0;
    mFragmentLength = 0;
    mReadGroup = null;
    mSuperCigar = null;
    mOverlapBases = null;
    mOverlapQuality = null;
    mOverlapInstructions = null;
    mCgReadDelta = null;
    mGenome = 0;
  }

  /**
   * Construct a new alignment record populated from a SAM record.
   * @param record SAM record. Requires header with sequence dictionary (for reference index lookup)
   * @param genome genome code for this record
   * @param chooser machine error chooser
   * @param minBaseQuality minimum read base quality
   */
  public VariantAlignmentRecord(final SAMRecord record, final int genome, MachineErrorChooserInterface chooser, int minBaseQuality) {
    super(record.getReferenceIndex(), record.getAlignmentStart() - 1, record.getReadUnmappedFlag() ? record.getAlignmentStart() - 1 + record.getReadLength() : record.getAlignmentEnd()); // picard end position is 1-based inclusive == 0-based exclusive
    mGenome = genome;
    mFragmentLength = record.getInferredInsertSize();
    final byte[] readBases = record.getReadBases();

    mBases = byteDNAtoByteHandleEquals(readBases); //we assume something will convert the = before use


    final byte[] baseQualities = record.getBaseQualities();
//    if (baseQualities.length == 0) {
//      mRecalibratedQuality = baseQualities;
//    } else {
//      mRecalibratedQuality = baseQualities;

    mCigar = record.getCigarString();
    mMappingQuality = (byte) record.getMappingQuality();
    mReadGroup = record.getReadGroup();
    mAmbiguity = MathUtils.unboxNatural(SamUtils.getNHOrIH(record));
    mAlignmentScore = MathUtils.unboxNatural(record.getIntegerAttribute("AS"));
    mSuperCigar = record.getStringAttribute(SamUtils.CG_SUPER_CIGAR);
    mMateSequenceId = record.getMateReferenceIndex();
    int f = 0;
    if (record.getReadPairedFlag()) {
      f += FLAG_PAIRED;
      if (record.getProperPairFlag()) {
        f += FLAG_MATED;
      }
      if (record.getFirstOfPairFlag()) {
        f += FLAG_FIRST;
      }
    }
    if (record.getReadNegativeStrandFlag()) {
      f += FLAG_NEGATIVE;
    }
    if (record.getReadUnmappedFlag()) {
      f += FLAG_UNMAPPED;
    }
    mFlag = (byte) f;

    final String overlapQuality = mSuperCigar == null
      ? SamUtils.allowEmpty(record.getStringAttribute(SamUtils.ATTRIBUTE_CG_OVERLAP_QUALITY))
      : SamUtils.allowEmpty(record.getStringAttribute(SamUtils.CG_SUPER_CIGAR_OVERLAP_QUALITY));
    final String cgOverlap = record.getStringAttribute(SamUtils.ATTRIBUTE_CG_OVERLAP_BASES);
    mOverlapBases = cgOverlap == null ? new byte[0] : cgOverlap.getBytes();
    DnaUtils.encodeArray(mOverlapBases);
    mOverlapInstructions = record.getStringAttribute(SamUtils.ATTRIBUTE_CG_RAW_READ_INSTRUCTIONS);
    mCgReadDelta = record.getStringAttribute(SamUtils.CG_READ_DELTA);

    // Perform recalibration at this point
    final PhredScaler me = chooser == null ? null : chooser.machineErrors(record.getReadGroup(), record.getReadPairedFlag());

    final int readLength = baseQualities.length + overlapQuality.length();

    final int backStepPosition;
    if (!overlapQuality.isEmpty()) {
      // Work out CG backstep position
      final boolean v1 = readLength == CgUtils.CG_RAW_READ_LENGTH;
      if (v1) {
        final boolean first = record.getReadPairedFlag() && record.getFirstOfPairFlag();
        backStepPosition = record.getReadNegativeStrandFlag() ^ !first ? baseQualities.length - CgUtils.CG_OVERLAP_POSITION : CgUtils.CG_OVERLAP_POSITION;
      } else {
        backStepPosition = record.getReadNegativeStrandFlag() ? baseQualities.length - CgUtils.CG2_OVERLAP_POSITION : CgUtils.CG2_OVERLAP_POSITION;
      }
    } else {
      backStepPosition = -1;
    }

    mRecalibratedQuality = new byte[baseQualities.length];
    final int machineStep;
    int machineCycle;
    if (record.getReadNegativeStrandFlag()) {
      machineCycle = readLength - 1;
      machineStep = -1;
    } else {
      machineCycle = 0;
      machineStep = 1;
    }
    int qualityPosition = 0;

    final Arm arm = !record.getReadPairedFlag() || record.getFirstOfPairFlag() ? Arm.LEFT : Arm.RIGHT;
    while (qualityPosition < backStepPosition && qualityPosition < mRecalibratedQuality.length) {
      final byte quality = baseQualities[qualityPosition];
      final int recalibrated = me == null ? quality : me.getScaledPhred(quality, machineCycle, arm);
      mRecalibratedQuality[qualityPosition] = (byte) recalibrated;
      machineCycle += machineStep;
      ++qualityPosition;
    }

    mOverlapQuality = new byte[overlapQuality.length()];
    for (int i = 0; qualityPosition < mRecalibratedQuality.length && i < overlapQuality.length(); ++i) {
      final byte scoreChar = (byte) (overlapQuality.charAt(i) - FastaUtils.PHRED_LOWER_LIMIT_CHAR);
      // Be careful to invoke the me.getScaledPhred that takes a char. It will correct for ascii encoding
      final int recalibrated = me == null ? scoreChar : me.getScaledPhred(scoreChar, machineCycle, arm);
      mOverlapQuality[i] = (byte) recalibrated;
      machineCycle += machineStep;
    }

    while (qualityPosition < mRecalibratedQuality.length) {
      final byte quality = baseQualities[qualityPosition];
      final int recalibrated = me == null ? quality : me.getScaledPhred(quality, machineCycle, arm);
      mRecalibratedQuality[qualityPosition] = (byte) recalibrated;
      machineCycle += machineStep;
      ++qualityPosition;
    }
    for (int i = 0; i < mRecalibratedQuality.length; ++i) {
      if (mRecalibratedQuality[i] < minBaseQuality) {
        if (MIN_QUALITY_AS_TWO) {
          mRecalibratedQuality[i] = 2;
        } else {
          mBases[i] = 0;
        }
      }
    }
  }

  /**
   * Test if this record represents an overflow condition.
   * @return if this is an overflow record
   */
  public boolean isOverflow() {
    return mFlag == -1;
  }

  /**
   * Construct a new alignment record populated from a SAM record.
   * @param record SAM record
   */
  public VariantAlignmentRecord(final SAMRecord record) {
    this(record, record.getReferenceIndex(), new DefaultMachineErrorChooser(), 0);
  }

  //XXX note we actually modify the contents of mBases in SuperCigarParser.updateReadWithTemplate
  public byte[] getRead() {
    return mBases;
  }

  /**
   * Get the binary phred quality values as a byte array (not Ascii).
   * @return quality
   */
  public byte[] getRecalibratedQuality() {
    return mRecalibratedQuality;
  }

  public String getCigar() {
    return mCigar;
  }

  @Override
  public int getMappingQuality() {
    return mMappingQuality & 0xFF;
  }

  // todo somehow make this an int
  public SAMReadGroupRecord getReadGroup() {
    return mReadGroup;
  }

  @Override
  public int getNHOrIH() {
    return mAmbiguity;
  }

  @Override
  public boolean isMated() {
    return (mFlag & FLAG_MATED) != 0;
  }

  public boolean isReadPaired() {
    return (mFlag & FLAG_PAIRED) != 0;
  }

  public boolean isFirst() {
    return (mFlag & FLAG_FIRST) != 0;
  }

  public Arm getArm() {
    return isFirst() ? Arm.LEFT : Arm.RIGHT;
  }

  public boolean isNegativeStrand() {
    return (mFlag & FLAG_NEGATIVE) != 0;
  }

  public boolean isUnmapped() {
    return (mFlag & FLAG_UNMAPPED) != 0;
  }

  // todo hopefully can get rid of this ... by rolling into normal cigar field
  public String getSuperCigar() {
    return mSuperCigar;
  }

  /** @return the quality of the overlapped bases as a byte[] (not ascii) */
  public byte[] getOverlapBases() {
    return mOverlapBases;
  }

  public byte[] getOverlapQuality() {
    return mOverlapQuality;
  }

  public String getOverlapInstructions() {
    return mOverlapInstructions;
  }

  public String getCGReadDelta() {
    return mCgReadDelta;
  }

  @Override
  public String toString() {
    return getStart() + " " + getCigar() + " " + DnaUtils.bytesToSequenceIncCG(getRead()) + " " + new String(FastaUtils.rawToAsciiQuality(getRecalibratedQuality()));
  }

  @Override
  public int disambiguateDuplicate(VariantAlignmentRecord rec) {
    final CompareHelper helper = new CompareHelper()
      .compare(getStart(), rec.getStart())
      .compare(getCigar(), rec.getCigar())
      .compare(compare(getRead(), rec.getRead()))
      .compare(compare(getRecalibratedQuality(), rec.getRecalibratedQuality()));
    return helper.result();
  }

  private static int compare(final byte[] a, final byte[] b) {
    if (a.length != b.length) {
      return a.length - b.length;
    }
    for (int k = 0; k < a.length; ++k) {
      if (a[k] != b[k]) {
        return a[k] - b[k];
      }
    }
    return 0;
  }

  private static int compare(final String a, final String b) {
    if (a != null) {
      return a.compareTo(b);
    } else if (b != null) {
      return -1;
    }
    return 0;
  }

  /**
   * Compares the object just on its values. the regular <code>compareTo</code> method will only return 0 for the same instance
   * @param var the object to compare
   * @return <code>-ve</code> for this object is less than, <code>0</code> for equal
   */
  public int valueCompareTo(VariantAlignmentRecord var) {
    // call starting at leftmost sequence, leftmost position, leftmost end
    if (var == this) {
      return 0;
    }
    final int thisRef = getSequenceId();
    final int thatRef = var.getSequenceId();
    if (thisRef == -1 && thatRef != -1) {
      return 1;
    } else if (thatRef == -1 && thisRef != -1) {
      return -1;
    }
    if (thisRef < thatRef) {
      return -1;
    }
    if (thisRef > thatRef) {
      return 1;
    }

    if (var.getStart() > getStart()) {
      return -1;
    } else if (var.getStart() < getStart()) {
      return 1;
    }

    final int aLen = getLength();
    final int bLen = var.getLength();
    if (aLen > bLen) {
      return 1;
    } else if (aLen < bLen) {
      return -1;
    }

    final int c = var.getCigar().compareTo(getCigar());
    if (c != 0) {
      return c;
    }
    final int mq = var.getMappingQuality() - getMappingQuality();
    if (mq != 0) {
      return mq;
    }
    final int r = compare(var.getRead(), getRead());
    if (r != 0) {
      return r;
    }
    final int q = compare(var.getRecalibratedQuality(), getRecalibratedQuality());
    if (q != 0) {
      return q;
    }
    final int f = var.mFlag - mFlag;
    if (f != 0) {
      return f;
    }
    final int as = var.mAlignmentScore - mAlignmentScore;
    if (as != 0) {
      return as;
    }
    final int nh = var.mAmbiguity - mAmbiguity;
    if (nh != 0) {
      return nh;
    }
    final int gg = var.mGenome - mGenome;
    if (gg != 0) {
      return gg;
    }
    final int sc = compare(var.mSuperCigar, mSuperCigar);
    if (sc != 0) {
      return sc;
    }
    return compare(Arrays.toString(var.mOverlapBases) + Arrays.toString(var.mOverlapQuality) + var.mOverlapInstructions + var.mCgReadDelta, Arrays.toString(mOverlapBases) + Arrays.toString(mOverlapQuality) + mOverlapInstructions + mCgReadDelta);
  }

  @Override
  public int compareTo(final VariantAlignmentRecord var) {
    final int ov = valueCompareTo(var);
    if (ov != 0) {
      return ov;
    }
    //    System.out.println(this + " cf. " + var);
    return System.identityHashCode(this) - System.identityHashCode(var);
  }

  @Override
  public boolean equals(final Object var) {
    return this == var;
  }

  @Override
  public int hashCode() {
    return super.hashCode();
  }

  /**
   * Genome this record corresponds to.
   * @return genome
   */
  @Override
  public int getGenome() {
    if (mGenome == -1) {
      throw new UnsupportedOperationException();
    }
    return mGenome;
  }

  private VariantAlignmentRecord mChainedRecord = null;

  @Override
  public void setNextInChain(final VariantAlignmentRecord rec) {
    mChainedRecord = rec;
  }

  @Override
  public VariantAlignmentRecord chain() {
    return mChainedRecord;
  }

  @Override
  public int getMateSequenceId() {
    return mMateSequenceId;
  }

  @Override
  public int getFragmentLength() {
    return mFragmentLength;
  }

  /**
   * Convert a byte array containing valid nucleotide characters into a byte array
   * with values (0=N, ... 4=T).
   * @param dna the string to be converted.
   * @return the byte array.
   */
  public static byte[] byteDNAtoByteHandleEquals(final byte[] dna) {
    final byte[] dnaBytes = new byte[dna.length];
    for (int i = 0; i < dna.length; ++i) {
      final char charAt = (char) dna[i];
      if (charAt == '=') {
        dnaBytes[i] = (byte) '=';
      } else {
        dnaBytes[i] = (byte) DNA.getDNA(charAt);
      }
    }
    return dnaBytes;
  }

}

