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
import java.util.HashMap;
import java.util.Map;

import com.rtg.ngs.PairedTopRandomImplementation.HitRecord;
import com.rtg.ngs.blocking.MapQScoringReadBlocker;
import com.rtg.ngs.blocking.ReadBlocker;
import com.rtg.ngs.tempstage.BinaryTempFileRecord;
import com.rtg.reader.PrereadNamesInterface;
import com.rtg.reader.SequencesReader;
import com.rtg.sam.SamBamConstants;
import com.rtg.sam.SamUtils;

import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMRecord;

/**
 * Concatenates and filters mated SAM records, filtering on the <code>XA</code> attribute.
 *
 */
public final class MatedSamResultsFilter extends AbstractSamResultsFilter {

  private final MapQScoringReadBlocker mXaBlocker;
  private final ReadBlocker mFreqBlockerLeft;
  private final ReadBlocker mFreqBlockerRight;
  private final SequencesReader mReader1;
  private final SequencesReader mReader2;

  private PairedTopRandomImplementation.HitRecord[] mHitsToKeep;
  private Map<String, Integer> mTemplateNames;

  /**
   * This version operates on mated results.
   *
   * @param xaBlocker blocker for combo score
   * @param freqBlockerLeft frequency blocker for left reads
   * @param freqBlockerRight frequency blocker for right reads
   * @param reader1 reader for left side
   * @param reader2 reader for right side (null if single end)
   * @param cg reads are from complete genomics
   * @param readIdOffset adjust read IDs by this much to put them in correct number space
   * @param readGroupId read group id
   * @param legacyCigars true if legacy cigar mode is enabled
   */
  public MatedSamResultsFilter(MapQScoringReadBlocker xaBlocker, ReadBlocker freqBlockerLeft, ReadBlocker freqBlockerRight, SequencesReader reader1, SequencesReader reader2, boolean cg, long readIdOffset, String readGroupId, boolean legacyCigars) {
    super(reader1, reader2, readGroupId, cg, readIdOffset, false, legacyCigars);
    mXaBlocker = xaBlocker;
    mFreqBlockerLeft = freqBlockerLeft;
    mFreqBlockerRight = freqBlockerRight;
    mReader1 = reader1;
    mReader2 = reader2;
    mHitsToKeep = null;
    mTemplateNames = null;
  }

  /* For toprandom */
  void setHitsToKeep(HitRecord[] hits) {
    mHitsToKeep = hits;
  }

  /* for toprandom */
  void setTemplateNames(PrereadNamesInterface names) {
    mTemplateNames = new HashMap<>();
    for (int i = 0; i < names.length(); i++) {
      mTemplateNames.put(names.name(i), i);
    }
  }

  @Override
  protected String getName() {
    return "Mated";
  }

  @Override
  protected SAMRecord filterRecord(SAMFileWriter samWriter, BinaryTempFileRecord rec, PrereadNamesInterface templateNames) throws IOException {
    assert rec.isReadPaired();
    final int readId = rec.getReadId();
    final int flag = rec.getSamFlags() & 0xff;
    assert (flag & SamBamConstants.SAM_READ_IS_PAIRED) != 0;
    assert (flag & SamBamConstants.SAM_READ_IS_MAPPED_IN_PROPER_PAIR) != 0;
    assert (flag & SamBamConstants.SAM_READ_IS_UNMAPPED) == 0;
    assert (flag & SamBamConstants.SAM_MATE_IS_UNMAPPED) == 0;
    final int score = rec.getComboScore();
    assert score >= 0;
    final boolean firstOfPair = (flag & SamBamConstants.SAM_READ_IS_FIRST_IN_PAIR) != 0;
    if (mHitsToKeep != null) {  // see if we output this record.
      final int templateId = rec.getReferenceId();
      final HitRecord curRec = mHitsToKeep[readId];
      if (curRec.mTemplateId != templateId) {
        return null;
      }
      final int templateStart;
      final int mateTemplateStart;
      final boolean reverse;
      final boolean mateReverse;
      if (firstOfPair) {
        templateStart = rec.getStartPosition();
        mateTemplateStart = rec.getMatePosition();
        reverse = (flag & SamBamConstants.SAM_READ_IS_REVERSE) != 0;
        mateReverse = (flag & SamBamConstants.SAM_MATE_IS_REVERSE) != 0;
      } else {
        templateStart = rec.getMatePosition();
        mateTemplateStart = rec.getStartPosition();
        reverse = (flag & SamBamConstants.SAM_MATE_IS_REVERSE) != 0;
        mateReverse = (flag & SamBamConstants.SAM_READ_IS_REVERSE) != 0;
      }
      if (templateStart - 1 != curRec.mZeroBasedTemplateStart || mateTemplateStart  - 1 != curRec.mZeroBasedMateTemplateStart || reverse != curRec.mReverse || mateReverse != curRec.mMateReverse) {
        return null;
      }
    }

    if (mXaBlocker.isBlocked1(readId, score) || mFreqBlockerLeft.isBlocked(readId) || mFreqBlockerRight.isBlocked(readId)) {
      //System.err.println("skipping ... " + readId);
      return null;
    } else {
      final SAMRecord record = new SAMRecord(mHeader);
      final boolean gq;
      final String readString;
      if (mCG) {
        readString = new String(rec.getCgReadString());
      } else {
        if (firstOfPair) {
          readString = read(readId, mReader1, (flag & SamBamConstants.SAM_READ_IS_REVERSE) != 0);
        } else {
          readString = read(readId, mReader2, (flag & SamBamConstants.SAM_READ_IS_REVERSE) != 0);
        }
      }
      if (firstOfPair) {
        assert (flag & SamBamConstants.SAM_READ_IS_SECOND_IN_PAIR) == 0;
        gq = qualField(readId, mReader1, (flag & SamBamConstants.SAM_READ_IS_REVERSE) != 0, rec, true, readString.length());
      } else {
        assert (flag & SamBamConstants.SAM_READ_IS_SECOND_IN_PAIR) != 0;
        gq = qualField(readId, mReader2, (flag & SamBamConstants.SAM_READ_IS_REVERSE) != 0, rec, false, readString.length());
      }
      final String readName = processReadName(mNames, readId, true);

      final int count = mXaBlocker.getCount1(readId);
      final int mapq = mXaBlocker.getMatedMapQ(readId);
      record.setReadName(readName);
      // Here we set the "alignment is secondary" flag if count > 1
      record.setFlags(count > 1 ? flag | SamBamConstants.SAM_SECONDARY_ALIGNMENT : flag);
      final int templateId = rec.getReferenceId();
      record.setReferenceName(templateNames.name(templateId));
      record.setAlignmentStart(rec.getStartPosition());
      record.setMappingQuality(mapq);
      record.setCigarString(new String(rec.getCigarString()));
      record.setReadString(readString);

      addBaseAttributes(record, rec);
      record.setMateReferenceIndex(record.getReferenceIndex());
      record.setMateAlignmentStart(rec.getMatePosition());
      record.setInferredInsertSize(rec.getTemplateLength());
      if (mQualBuffer != null) {
        final byte[] quals = new byte[mQualStringLength];
        System.arraycopy(mFlattenedQualBuffer, 0, quals, 0, mQualStringLength);
        record.setBaseQualities(quals);
        if (gq) {
          final String cgQual = new String(mGQBuffer, 0, mGQLength);
          record.setAttribute(SamUtils.CG_SUPER_CIGAR_OVERLAP_QUALITY, cgQual);
        }
      }
      assert count > 0 : mXaBlocker.isBlocked1(readId, score) + " " + mFreqBlockerLeft.isBlocked(readId) + " " + mFreqBlockerRight.isBlocked(readId);
      addIhAttributes(record, readId, mHitsToKeep == null, count, mXaBlocker);
      samWriter.addAlignment(record);

      return record;
    }
  }
}
