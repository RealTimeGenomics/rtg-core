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

import com.rtg.ngs.SingleEndTopRandomImplementation.HitRecord;
import com.rtg.ngs.blocking.MapQScoringReadBlocker;
import com.rtg.ngs.blocking.ReadBlocker;
import com.rtg.ngs.tempstage.BinaryTempFileRecord;
import com.rtg.pairedend.ReadStatusListener;
import com.rtg.reader.NamesInterface;
import com.rtg.reader.PrereadType;
import com.rtg.reader.SequencesReader;
import com.rtg.sam.SamBamConstants;

import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMRecord;

/**
 * Concatenates one or more intermediate single end SAM files and filters out just the
 * records with the best score. It also adds the <code>HI</code> and <code>NH/IH</code>
 * attributes to each output record.
 *
 */
public final class SingleEndSamResultsFilter extends AbstractSamResultsFilter {

  private final MapQScoringReadBlocker mAsBlocker;
  private final ReadBlocker mFreqBlocker;
  private final ReadStatusListener mListener;
  private HitRecord[] mHitsToKeep;
  private final SequencesReader mReader;

  /**
   * This version operates on single ended results, it operates on the AS field.
   *
   * @param asBlocker information about the min score for each read.
   * @param freqBlocker information about the repeat frequence for each read.
   * @param listener each record written will set a flag in this listener
   * @param readIdOffset adjust read IDs by this much to put them in the correct range
   * @param reader1 reader for left side
   * @param readGroupId the id of the read group for this run (may be null)
   * @param legacyCigars true if legacy cigar mode is enabled
   */
  public SingleEndSamResultsFilter(MapQScoringReadBlocker asBlocker, ReadBlocker freqBlocker, ReadStatusListener listener, long readIdOffset, SequencesReader reader1, String readGroupId, boolean legacyCigars) {
    super(reader1, null, readGroupId, reader1 != null && PrereadType.CG == reader1.getPrereadType(), readIdOffset, false, legacyCigars);
    mAsBlocker = asBlocker;
    mFreqBlocker = freqBlocker;
    mListener = listener;
    mReader = reader1;
  }

  @Override
  protected String getName() {
    return "Alignment";
  }

  /* For toprandom */
  void setHitsToKeep(HitRecord[] hits) {
    mHitsToKeep = hits;
  }

  /* for toprandom */
  void setTemplateNames(NamesInterface names) {
  }

  @Override
  protected SAMRecord filterRecord(SAMFileWriter samWriter, BinaryTempFileRecord rec, NamesInterface templateNames) throws IOException {
    final int readId = rec.getReadId();
    final int flag = rec.getSamFlags() & 0xff;
    assert (flag & SamBamConstants.SAM_READ_IS_PAIRED) == 0;
    assert (flag & SamBamConstants.SAM_READ_IS_UNMAPPED) == 0;
    final int score = rec.getAlignmentScore();
    assert score >= 0;

    if (mHitsToKeep != null) {
      final HitRecord curRec = mHitsToKeep[readId];
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
//    System.err.println("readId=" + readId + ", score=" + score);
    if (!mAsBlocker.isBlocked1(readId, score) && !mFreqBlocker.isBlocked(readId)) {
      final SAMRecord record = new SAMRecord(mHeader);

      final int count = mAsBlocker.getCount1(readId);
      final int mapq = mAsBlocker.getMapQ(readId);
      record.setReadName(processReadName(mNames, readId, false));
      // Here we set the "alignment is secondary" flag if count > 1
      record.setFlags(count > 1 ? flag | SamBamConstants.SAM_SECONDARY_ALIGNMENT : flag);

      record.setReferenceName(templateNames.name(rec.getReferenceId()));

      record.setAlignmentStart(rec.getStartPosition());
      record.setMappingQuality(mapq);
      record.setCigarString(new String(rec.getCigarString()));
      final String readString;
      readString = read(readId, mReader, (flag & SamBamConstants.SAM_READ_IS_REVERSE) != 0);
      qualField(readId, mReader, (flag & SamBamConstants.SAM_READ_IS_REVERSE) != 0, rec, true, readString.length());
      if (mQualBuffer != null) {
        final byte[] quals = new byte[mQualStringLength];
        System.arraycopy(mFlattenedQualBuffer, 0, quals, 0, mQualStringLength);
        record.setBaseQualities(quals);
      }
      record.setReadString(readString);

      addBaseAttributes(record, rec);
      addIhAttributes(record, readId, mHitsToKeep == null, count, mAsBlocker);

      samWriter.addAlignment(record);

      if (count == 1) {
        mListener.addStatus(readId, ReadStatusTracker.UNIQUELY_MAPPED_FIRST);
      }
      mListener.addStatus(readId, ReadStatusTracker.UNMATED_FIRST);
      return record;
    }
    return null;
  }
}
