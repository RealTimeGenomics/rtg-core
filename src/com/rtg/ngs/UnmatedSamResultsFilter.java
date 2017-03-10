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
package com.rtg.ngs;

import java.io.IOException;

import com.rtg.index.hash.ngs.ReadEncoder;
import com.rtg.ngs.SingleEndTopRandomImplementation.HitRecord;
import com.rtg.ngs.blocking.MapQScoringReadBlocker;
import com.rtg.ngs.blocking.ReadBlocker;
import com.rtg.ngs.tempstage.BinaryTempFileRecord;
import com.rtg.reader.NamesInterface;
import com.rtg.reader.PrereadType;
import com.rtg.reader.SequencesReader;
import com.rtg.sam.SamBamConstants;
import com.rtg.sam.SamUtils;
import com.rtg.variant.sv.UnmatedAugmenter;

import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMRecord;

/**
 * Concatenates and filters unmated SAM files, filtering on the <code>AS</code> attribute.
 *
 */
public final class UnmatedSamResultsFilter extends AbstractSamResultsFilter {

  private final MapQScoringReadBlocker mAsBlockerLeft;
  private final MapQScoringReadBlocker mAsBlockerRight;
  private final ReadBlocker mFreqBlockerLeft;
  private final ReadBlocker mFreqBlockerRight;
  private HitRecord[] mHitsToKeep = null;

  private final UnmatedAugmenter mAugmenter;

  private SequencesReader mReader1;
  private SequencesReader mReader2;
  private boolean mCG;

  /**
   * This version operates on unmated results.
   *
   * @param asBlockerLeft information about the min score for each read for the left hand side.
   * @param asBlockerRight information about the min score for each read for the right hand side.
   * @param freqBlockerLeft frequency blocker for left reads
   * @param freqBlockerRight frequency blocker for right reads
   * @param readIdOffset adjust the read IDs by this much to put them in the correct range
   * @param reader1 reader for left side
   * @param reader2 reader for right side (null if single end)
   * @param readGroupId the id of the read group for this run (may be null)
   * @param legacyCigars true if legacy cigar mode is enabled
   * @param augmenterMerger supplier of unmated augmenters
   */
  public UnmatedSamResultsFilter(MapQScoringReadBlocker asBlockerLeft, MapQScoringReadBlocker asBlockerRight, ReadBlocker freqBlockerLeft, ReadBlocker freqBlockerRight, long readIdOffset, SequencesReader reader1, SequencesReader reader2, String readGroupId, boolean legacyCigars, UnmatedAugmenter.Merger augmenterMerger) {
    super(reader1, reader2, readGroupId, PrereadType.CG.equals(reader1.getPrereadType()), readIdOffset, false, legacyCigars);
    mAsBlockerLeft = asBlockerLeft;
    mAsBlockerRight = asBlockerRight;
    mFreqBlockerLeft = freqBlockerLeft;
    mFreqBlockerRight = freqBlockerRight;
    mReader1 = reader1;
    mReader2 = reader2;
    mAugmenter = augmenterMerger == null ? null : augmenterMerger.createUnmatedAugmenter();

  }

  /* For toprandom */
  void setHitsToKeep(HitRecord[] hits) {
    mHitsToKeep = hits;
  }

  /* for toprandom */
  void setTemplateNames(NamesInterface names) {
  }

  @Override
  protected String getName() {
    return "Unmated";
  }


  @Override
  protected SAMRecord filterRecord(SAMFileWriter samWriter, BinaryTempFileRecord rec, NamesInterface templateNames) throws IOException {
    final int readId = rec.getReadId();
    final int flag = rec.getSamFlags() & 0xff;
    assert (flag & SamBamConstants.SAM_READ_IS_PAIRED) != 0;
    assert (flag & SamBamConstants.SAM_READ_IS_MAPPED_IN_PROPER_PAIR) == 0;
    assert (flag & SamBamConstants.SAM_READ_IS_UNMAPPED) == 0;
    assert (flag & SamBamConstants.SAM_MATE_IS_UNMAPPED) == 0;  // We will set it in this class
    assert (flag & SamBamConstants.SAM_READ_IS_FIRST_IN_PAIR) != 0 || (flag & SamBamConstants.SAM_READ_IS_SECOND_IN_PAIR) != 0;
    final int score = rec.getAlignmentScore();
    assert score >= 0;

    if (mHitsToKeep != null) {
      final int encodedReadId;
      if ((flag & SamBamConstants.SAM_READ_IS_FIRST_IN_PAIR) != 0) {
        encodedReadId = ReadEncoder.PAIRED_FIRST.encode(readId);
      } else {
        encodedReadId = ReadEncoder.PAIRED_SECOND.encode(readId);
      }

      final HitRecord curRec = mHitsToKeep[encodedReadId];
      final int templateId = rec.getReferenceId();
      if (curRec.mTemplateId != templateId) {
        return null;
      }
      final int templateStart = rec.getStartPosition();
      final boolean reverse = (flag & SamBamConstants.SAM_READ_IS_REVERSE) != 0;
      if (templateStart - 1 != curRec.mZeroBasedTemplateStart
          || reverse != curRec.mReverse) {
        return null;
      }
    }
    assert mPaired;
    //System.out.println("readId=" + readId + ", score=" + score);
    if ((flag & SamBamConstants.SAM_READ_IS_FIRST_IN_PAIR) != 0) {
      assert (flag & SamBamConstants.SAM_READ_IS_SECOND_IN_PAIR) == 0;
    } else {
      assert (flag & SamBamConstants.SAM_READ_IS_SECOND_IN_PAIR) != 0;
    }
    final String readName = processReadName(mNames, readId, true);

    final boolean first = (flag & SamBamConstants.SAM_READ_IS_FIRST_IN_PAIR) != 0;
    final MapQScoringReadBlocker asblocker = first ? mAsBlockerLeft : mAsBlockerRight;
    final ReadBlocker freqblocker = first ? mFreqBlockerLeft : mFreqBlockerRight;
    final boolean keepRecord = !asblocker.isBlocked1(readId, score) && !freqblocker.isBlocked(readId);
    if (keepRecord) {
      final SAMRecord record = new SAMRecord(mHeader);
      final int count = asblocker.getCount1(readId);
      final int mapq = asblocker.getMapQ(readId);
      record.setReadName(readName);

      // Here we set the "alignment is secondary" flag if count > 1
      int newflag = flag;
      final MapQScoringReadBlocker mateasblocker = !first ? mAsBlockerLeft : mAsBlockerRight;
      final ReadBlocker matefreqblocker = !first ? mFreqBlockerLeft : mFreqBlockerRight;
      final boolean isMateMapped = (mateasblocker.getCount1(readId) > 0) && !mateasblocker.isBlocked1(readId, mateasblocker.getScore1(readId)) && !matefreqblocker.isBlocked(readId);
      if (!isMateMapped) {
        newflag |= SamBamConstants.SAM_MATE_IS_UNMAPPED;
      }
      if (count > 1) {
        newflag |= SamBamConstants.SAM_SECONDARY_ALIGNMENT;
      }
      record.setFlags(newflag);

      final int templateId = rec.getReferenceId();
      record.setReferenceName(templateNames.name(templateId));

      record.setAlignmentStart(rec.getStartPosition());
      record.setMappingQuality(mapq);

      record.setCigarString(new String(rec.getCigarString()));
      final String readString;
      final boolean gq;
      final boolean firstOfPair = (flag & SamBamConstants.SAM_READ_IS_FIRST_IN_PAIR) != 0;
      if (mCG) { //if CG or no reads specified then read should be in the original SAM file
        readString = new String(rec.getCgReadString());
      } else if (firstOfPair) {
        readString = read(readId, mReader1, (flag & SamBamConstants.SAM_READ_IS_REVERSE) != 0);
      } else {
        readString = read(readId, mReader2, (flag & SamBamConstants.SAM_READ_IS_REVERSE) != 0);
      }
      if (firstOfPair) {
        gq = qualField(readId, mReader1, (flag & SamBamConstants.SAM_READ_IS_REVERSE) != 0, rec, true, readString.length());
      } else {
        gq = qualField(readId, mReader2, (flag & SamBamConstants.SAM_READ_IS_REVERSE) != 0, rec, false, readString.length());
      }
      addBaseAttributes(record, rec);
      if (mQualBuffer != null) {
        record.setMateAlignmentStart(rec.getMatePosition());
        record.setInferredInsertSize(rec.getTemplateLength());
        final byte[] quals = new byte[mQualStringLength];
        System.arraycopy(mFlattenedQualBuffer, 0, quals, 0, mQualStringLength);
        record.setBaseQualities(quals);
        if (gq) {
          final String cgQual = new String(mGQBuffer, 0, mGQLength);
          record.setAttribute(SamUtils.CG_SUPER_CIGAR_OVERLAP_QUALITY, cgQual);
        }
      }
      record.setReadString(readString);

      record.setMateReferenceIndex(-1);
      addIhAttributes(record, readId, mHitsToKeep == null, count, asblocker);

      samWriter.addAlignment(record);

      // add to unmated augmenter
      if (mAugmenter != null) {
        mAugmenter.processRecord(record);
      }

      return record;
    }
    return null;
  }

  void setReadersAndType(SequencesReader copy, SequencesReader copy2, boolean cg) {
    mReader1 = copy;
    mReader2 = copy2;
    final int maxLen = getQualLength(copy, copy2);
    mReadBuffer = new byte[maxLen];
    mCG = cg;
  }
}
