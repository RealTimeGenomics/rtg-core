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
package com.rtg.variant.cnv;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;

import com.rtg.launcher.NoStatistics;
import com.rtg.launcher.ParamsTask;
import com.rtg.reader.SequencesReader;
import com.rtg.sam.RecordIterator;
import com.rtg.sam.SamReadingContext;
import com.rtg.sam.SamRecordPopulator;
import com.rtg.sam.SamUtils;
import com.rtg.sam.ThreadedMultifileIterator;
import com.rtg.util.SingletonPopulatorFactory;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.ErrorType;
import com.rtg.util.diagnostic.SlimException;
import com.rtg.util.intervals.RegionRestriction;
import com.rtg.variant.cnv.region.Region;
import com.rtg.variant.cnv.region.RegionUtils;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;

/**
 */
public class CnvProductTask extends ParamsTask<CnvProductParams, NoStatistics> {

  private int[] mSequenceLengths;

  /**
   * Constructor
   *
   * @param params the params
   * @param reportStream the stream for reporting things like statistics
   */
  public CnvProductTask(CnvProductParams params, OutputStream reportStream) {
    super(params, reportStream, new NoStatistics(), null);
    mSequenceLengths = null;
    if (mParams.filterParams().restriction() != null && mParams.filterParams().restriction().getStart() != RegionRestriction.MISSING) {
      Diagnostic.warning("WARNING: using only a sub-region of the reads in a template sequence will not produce optimal results.");
    }
  }

  @Override
  protected void exec() throws IOException {
    final ArrayList<File> allFiles = new ArrayList<>(mParams.mappedBase());
    allFiles.addAll(mParams.mappedTarget());
    final SequencesReader reference = (mParams.genome() == null) ? null : mParams.genome().reader();
    final SAMFileHeader header = SamUtils.getUberHeader(reference, allFiles, mParams.ignoreIncompatibleSamHeaders(), null);
    final SAMSequenceDictionary dict = header.getSequenceDictionary();
    final Map<Integer, String> templateNameMap = makeTemplateNameMap(dict);
    setSequenceLengths(makeTemplateLengths(dict));
    final int[][] chunksBase;
    final SingletonPopulatorFactory<SAMRecord> pf = new SingletonPopulatorFactory<>(new SamRecordPopulator());
    final SamReadingContext cBase = new SamReadingContext(mParams.mappedBase(), mParams.threads(), mParams.filterParams(), header, reference);
    try (RecordIterator<SAMRecord> multiIteratorBase = new ThreadedMultifileIterator<>(cBase, pf)) {
      chunksBase = chunkSamFile(multiIteratorBase);
    }
    final int[][] chunksTarget;
    final SamReadingContext cTarget = new SamReadingContext(mParams.mappedTarget(), mParams.threads(), mParams.filterParams(), header, reference);
    try (RecordIterator<SAMRecord> multiIteratorTarget = new ThreadedMultifileIterator<>(cTarget, pf)) {
      chunksTarget = chunkSamFile(multiIteratorTarget);
    }
    final Map<String, Region> nregions;
    if (reference == null) {
      nregions = new HashMap<>();
    } else {
      nregions = RegionUtils.regionsFromSDF(reference, mParams.bucketSize());
    }
    new CnvRatio(mParams.magicConstant(), nregions, templateNameMap, mParams, (double) mSumLengths / mRecordCount, mParams.extraPenaltyOff()).exec(chunksBase, chunksTarget);
  }

  void setSequenceLengths(int[] lengths) {
    mSequenceLengths = lengths;
  }

  long mSumLengths = 0;
  long mRecordCount = 0;

  int[][] chunkSamFile(RecordIterator<SAMRecord> iterator) {
    final int[][] ret = new int[mSequenceLengths.length][];
    int[] currentStartCoverage = null;
    int currentSequenceId = -1;
    while (iterator.hasNext()) {
      final SAMRecord rec = iterator.next();
      if (rec.getReferenceIndex() != currentSequenceId) {
        fillBucket(currentStartCoverage, ret, currentSequenceId);
        currentSequenceId = rec.getReferenceIndex();
        final int bucketSize = Math.max((mSequenceLengths[currentSequenceId] + mParams.bucketSize() - 1) / mParams.bucketSize(), 1);
        ret[currentSequenceId] = new int[bucketSize];
        final int length = mSequenceLengths[currentSequenceId];
        if (length < 0) {
          throw new SlimException(ErrorType.INFO_ERROR, "Sequence length " + length + " is negative for sequence ID " + currentSequenceId);
        }
        currentStartCoverage = new int[length];
      }
      final int pos = rec.getAlignmentStart() - 1;
      if (pos < 0 || pos >= currentStartCoverage.length) {
        throw new SlimException(ErrorType.INFO_ERROR, "SAM Reading Error, pos: " + pos + ", sequenceIndex: " + rec.getReferenceIndex() + ", sequenceName: " + rec.getReferenceName() + ", coverage array size: " + currentStartCoverage.length + ", reflected array size: " + java.lang.reflect.Array.getLength(currentStartCoverage));
      } else {
        currentStartCoverage[pos]++;
        mSumLengths += rec.getReadLength();
        ++mRecordCount;
      }
    }
    fillBucket(currentStartCoverage, ret, currentSequenceId);
    return ret;
  }

  void fillBucket(int[] startCoverage, int[][] buckets, int sequenceId) {
    if (startCoverage == null || sequenceId < 0) {
      return;
    }
    int bucketIndex = 0;
    for (int i = 0, j = 0; i < startCoverage.length; ++i, ++j) {
      if (j == mParams.bucketSize()) {
        j = 0;
        ++bucketIndex;
      }
      if (startCoverage[i] > 0) {
        buckets[sequenceId][bucketIndex] += mParams.filterStartPositions() ? 1 : startCoverage[i];
      }
    }
  }

  static Map<Integer, String> makeTemplateNameMap(final SAMSequenceDictionary dictionary) {
    final Map<Integer, String> map = new HashMap<>(dictionary.size());
    for (final SAMSequenceRecord rec : dictionary.getSequences()) {
      map.put(rec.getSequenceIndex(), rec.getSequenceName());
    }
    return map;
  }

  static int[] makeTemplateLengths(final SAMSequenceDictionary dictionary) {
    final int[] lengths = new int[dictionary.size()];
    for (final SAMSequenceRecord rec : dictionary.getSequences()) {
      lengths[rec.getSequenceIndex()] = rec.getSequenceLength();
    }
    return lengths;
  }

}
