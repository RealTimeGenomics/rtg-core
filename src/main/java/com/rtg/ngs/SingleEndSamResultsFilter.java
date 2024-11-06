/*
 * Copyright (c) 2018. Real Time Genomics Limited.
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the
 *    distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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
