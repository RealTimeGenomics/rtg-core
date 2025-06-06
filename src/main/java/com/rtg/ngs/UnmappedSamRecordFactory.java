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

import java.io.Closeable;
import java.io.IOException;

import com.reeltwo.jumble.annotations.TestClass;
import com.rtg.mode.DnaUtils;
import com.rtg.reader.NamesInterface;
import com.rtg.reader.ReadHelper;
import com.rtg.reader.SequencesReader;
import com.rtg.reference.ReferenceGenome;
import com.rtg.sam.ReadGroupUtils;
import com.rtg.sam.SamUtils;
import com.rtg.util.intervals.LongRange;
import com.rtg.variant.sv.ReadGroupStatsCalculator;
import com.rtg.variant.sv.UnmatedAugmenter;

import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;

/**
 */
@TestClass("com.rtg.ngs.UnmappedSamAlignmentWriterTest")
public class UnmappedSamRecordFactory implements Closeable {

  private final SharedResources mSharedResources;
  SequencesReader mReader1;
  SequencesReader mReader2;
  private final NamesInterface mNames;
  private final LongRange mReadRange;
  private final SAMReadGroupRecord mReadGroupRecord;

  private final boolean mPaired;
  private UnmatedAugmenter mAugmenter = null;
  private ReadGroupStatsCalculator mCalculator = null;
  private MapReportData mMapReportData;
  private final ReferenceGenome mReferenceGenome;

    /**
   * Construct an unmapped record writer that can include unmated augmentation information
   *
   * @param params various parameters
   * @param sharedResources the shared resources
   * @throws IOException if an IO error occurs when reading preread names
   */
  public UnmappedSamRecordFactory(NgsParams params, SharedResources sharedResources) throws IOException {

    mSharedResources = sharedResources;
    mReader1 = mSharedResources.firstReaderCopy();
    mReader2 = mSharedResources.secondReaderCopy();

    if (params.outputParams().outputReadNames()) {
      mNames = mReader1.names();
    } else {
      mNames = null;
    }
    mReadRange = params.buildFirstParams().readerRestriction();
    mReadGroupRecord = params.outputParams().readGroup();
    mPaired = params.paired();
    final SequencesReader referenceReader = params.searchParams().reader();
    mReferenceGenome = ReferenceGenome.hasReferenceFile(referenceReader) ? new ReferenceGenome(referenceReader, params.sex()) : null;
  }

  /**
   * @param paired true if paired
   * @param readRange range containing read id offset
   * @param names source of names
   * @param readId internal read id
   * @return the read name
   */
  static String extractReadName(boolean paired, LongRange readRange, NamesInterface names, int readId) {
    if (names == null) {
      final long offset = Math.max(readRange.getStart(), 0);
      return String.valueOf(readId + offset);
    }
    return SamUtils.samReadName(names.name(readId), paired);
  }

  // If set, mapping statistics will include unmapped counts.
  void setMapReportData(MapReportData.Merger merger) {
    mMapReportData = merger == null ? null : merger.createMapReportData();
  }


  /**
   * Set information needed for augmenting unmapped records with unmated information.
   * @param augmenterMerger the unmated augmenter that contains all known pairing information for adding to unmapped reads
   * @param statsMerger the read group statistics calculator containing calculated insert sizes for use with the augmenter
   */
  void setAugmenterInfo(UnmatedAugmenter.Merger augmenterMerger, ReadGroupStatsCalculator.Merger statsMerger) {
    if (augmenterMerger == null) {
      mAugmenter = null;
    } else {
      mAugmenter = augmenterMerger.blend();
      mAugmenter.addMachineType(mReadGroupRecord);
    }
    mCalculator = statsMerger == null ? null : statsMerger.blend();
  }


  /**
   * Creates an unmapped SAM record
   * @param readId id of read
   * @param first true for first of pair.
   * @param xc value for <code>XC</code> attribute, <code>'\0'</code> for none
   * @param mateUnmapped true if the mate has not been mapped
   * @return the created SAM record
   */
  public SAMRecord unmappedResult(int readId, boolean first, char xc, boolean mateUnmapped) {
    final SequencesReader reader = !mPaired || first ? mReader1 : mReader2;
    final byte[] read = ReadHelper.getRead(reader, readId);
    final byte[] qual = ReadHelper.getQual(reader, readId);
    final SAMRecord rec = createUnmappedSamRecord(readId, DnaUtils.bytesToSequenceIncCG(read), qual, mPaired, !mPaired || first, xc, mPaired && mateUnmapped);
    if (mMapReportData != null) {
      mMapReportData.processRead(rec);
    }
    return rec;
  }


  protected SAMRecord createUnmappedSamRecord(int readId, String read, byte[] qual, boolean isPaired, boolean isFirst, char xc, boolean mateUnmapped) {
    final SAMRecord ret = new SAMRecord(mSharedResources.getHeader());
    final String readName;
    if (mNames != null) {
      readName = extractReadName(mPaired, mReadRange, mNames, readId);
    } else {
      readName = extractReadName(mPaired, mReadRange, null, readId);
    }
    ret.setReadName(readName);
    ret.setReferenceName("*");
    ret.setReadUnmappedFlag(true);
    ret.setReadPairedFlag(isPaired);
    ret.setReadBases(read.getBytes());
    if (qual != null) {
      ret.setBaseQualities(qual);
    }
    ret.setInferredInsertSize(0);
    ret.setMateAlignmentStart(0);
    ret.setMateReferenceName("*");
    if (isPaired) {
      ret.setFirstOfPairFlag(isFirst);
      ret.setSecondOfPairFlag(!isFirst);
      ret.setMateUnmappedFlag(mateUnmapped);
    }
    if ('\0' != xc) {
      ret.setAttribute(SamUtils.ATTRIBUTE_UNMAPPED_ETYMOLOGY, xc);
    }
    addSamRg(ret);

    if (mAugmenter != null) {
      mAugmenter.updateUnmappedRecord(ret, mCalculator, mReferenceGenome);
    }
    return ret;
  }

  private void addSamRg(final SAMRecord rec) {
    if (mReadGroupRecord != null) {
      rec.setAttribute(ReadGroupUtils.RG_ATTRIBUTE, mReadGroupRecord.getReadGroupId());
    }
  }

  @Override
  public void close() throws IOException {
    if (mReader1 != null) {
      mReader1.close();
      mReader1 = null;
    }
    if (mReader2 != null) {
      mReader2.close();
      mReader2 = null;
    }
  }
}
