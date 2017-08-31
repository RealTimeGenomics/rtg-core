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
