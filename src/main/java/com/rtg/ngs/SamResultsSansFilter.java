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

import com.rtg.ngs.tempstage.BinaryTempFileRecord;
import com.rtg.pairedend.ReadStatusListener;
import com.rtg.reader.NamesInterface;
import com.rtg.reader.PrereadType;
import com.rtg.reader.SequencesReader;
import com.rtg.sam.SamBamConstants;
import com.rtg.sam.SamUtils;

import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMRecord;

/**
 * Concatenates and filters unmated SAM files, filtering on absolutely nothing.
 */
public final class SamResultsSansFilter extends AbstractSamResultsFilter {

  private final ReadStatusListener mListener;
  private final SequencesReader mReader1;
  private final SequencesReader mReader2;


  /**
   * This version operates on whatever results.
   * @param listener each record written will set a flag in this listener
   * @param readIdOffset read offset
   * @param reader1 reader for left side
   * @param reader2 reader for right side (null if single end)
   * @param readGroupId the id of the read group for this run (may be null)
   * @param legacyCigars true if legacy cigar mode is enabled
   */
  public SamResultsSansFilter(ReadStatusListener listener, long readIdOffset, SequencesReader reader1, SequencesReader reader2, String readGroupId, boolean legacyCigars) {
    super(reader1, reader2, readGroupId, PrereadType.CG == reader1.getPrereadType(), readIdOffset, true, legacyCigars);
    mListener = listener;
    mReader1 = reader1;
    mReader2 = reader2;
  }

  @Override
  protected String getName() {
    return "Sans";
  }

  @Override
  protected SAMRecord filterRecord(SAMFileWriter samWriter, BinaryTempFileRecord rec, NamesInterface templateNames) throws IOException {
    final int readId = rec.getReadId();
    final int flag = rec.getSamFlags() & 0xff;
    //assert (flag & SamBamConstants.SAM_READ_IS_PAIRED) != 0;

    final boolean isPaired = (flag & SamBamConstants.SAM_READ_IS_PAIRED) != 0;
    assert isPaired == mPaired;
    assert !isPaired || (flag & SamBamConstants.SAM_READ_IS_MAPPED_IN_PROPER_PAIR) == 0;
    assert (flag & SamBamConstants.SAM_READ_IS_UNMAPPED) == 0;
    assert !isPaired || (flag & SamBamConstants.SAM_MATE_IS_UNMAPPED) == 0;
    assert !isPaired || (flag & SamBamConstants.SAM_READ_IS_FIRST_IN_PAIR) != 0 || (flag & SamBamConstants.SAM_READ_IS_SECOND_IN_PAIR) != 0;
    final int score = rec.getAlignmentScore();
    if (score < 0) {
      // oops! no AS:score in this record
      throw new RuntimeException("SAM record has no " + SamUtils.ATTRIBUTE_ALIGNMENT_SCORE + ":score field!  readId=" + readId + " flag=" + flag);
    }
    //System.out.println("readId=" + readId + ", score=" + score);
    final boolean firstOfPair = (flag & SamBamConstants.SAM_READ_IS_FIRST_IN_PAIR) != 0;
    if ((flag & SamBamConstants.SAM_READ_IS_PAIRED) == 0 || firstOfPair) {
      assert (flag & SamBamConstants.SAM_READ_IS_SECOND_IN_PAIR) == 0;
    } else {
      assert (flag & SamBamConstants.SAM_READ_IS_SECOND_IN_PAIR) != 0;
    }
    final String readName = processReadName(mNames, readId, mPaired);
    final SAMRecord record = new SAMRecord(mHeader);
    record.setReadName(readName);
    record.setFlags(flag);
    record.setAlignmentStart(rec.getStartPosition());

    // we don't do anything in regards to alignment is secondary flag
    final int templateId = rec.getReferenceId();
    record.setReferenceName(templateNames.name(templateId));
    record.setMappingQuality(255);
    record.setCigarString(new String(rec.getCigarString()));
    final boolean mated = (rec.getSamFlags() & SamBamConstants.SAM_READ_IS_MAPPED_IN_PROPER_PAIR) != 0;
    if (mated) {
      record.setMateReferenceIndex(record.getReferenceIndex());
      record.setMateAlignmentStart(rec.getMatePosition());
      record.setInferredInsertSize(rec.getTemplateLength());
    } else {
      record.setMateReferenceIndex(-1);
    }
    final String readString;
    if (mCG) {
      readString = new String(rec.getCgReadString());
    } else {
      if (!isPaired || firstOfPair) {
        readString = read(readId, mReader1, (flag & SamBamConstants.SAM_READ_IS_REVERSE) != 0);
      } else {
        readString = read(readId, mReader2, (flag & SamBamConstants.SAM_READ_IS_REVERSE) != 0);
      }
    }
    final boolean gq;
    if (!isPaired || firstOfPair) {
      gq = qualField(readId, mReader1, (flag & SamBamConstants.SAM_READ_IS_REVERSE) != 0, rec, true, readString.length());
    } else {
      gq = qualField(readId, mReader2, (flag & SamBamConstants.SAM_READ_IS_REVERSE) != 0, rec, false, readString.length());
    }
    if (mQualBuffer != null) {
      final byte[] quals = new byte[mQualStringLength];
      System.arraycopy(mFlattenedQualBuffer, 0, quals, 0, mQualStringLength);
      record.setBaseQualities(quals);
      if (gq) {
        final String cgQual = new String(mGQBuffer, 0, mGQLength);
        record.setAttribute(SamUtils.CG_SUPER_CIGAR_OVERLAP_QUALITY, cgQual);
      }
    }
    record.setReadString(readString);
    addBaseAttributes(record, rec);
    if (mPaired) {
      if (rec.isUnfilteredMated()) {
        record.setAttribute("XM", 1);
      }
    }
    samWriter.addAlignment(record);

    if ((flag & SamBamConstants.SAM_READ_IS_PAIRED) == 0 || (flag & SamBamConstants.SAM_READ_IS_FIRST_IN_PAIR) != 0) {
      mListener.addStatus(readId, ReadStatusTracker.UNMATED_FIRST);
    } else {
      mListener.addStatus(readId, ReadStatusTracker.UNMATED_SECOND);
    }
    return record;
  }
}
