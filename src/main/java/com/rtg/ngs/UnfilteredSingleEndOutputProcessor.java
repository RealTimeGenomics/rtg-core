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

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.util.ArrayList;

import com.rtg.index.hash.ngs.OutputProcessor;
import com.rtg.launcher.HashingRegion;
import com.rtg.ngs.blocking.MapQScoringReadBlocker;
import com.rtg.ngs.blocking.ReadBlocker;
import com.rtg.ngs.tempstage.BinaryTempFileRecord;
import com.rtg.ngs.tempstage.SingleEndTempFileWriter;
import com.rtg.reader.NamesInterface;
import com.rtg.util.NullStreamUtils;
import com.rtg.util.Pair;
import com.rtg.util.StringUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.io.FileUtils;

/**
 * Quick and dirty writer of unfiltered SAM output files.
 *
 */
public class UnfilteredSingleEndOutputProcessor extends AbstractSdfOutputProcessor {

  private static final int MATCHED = ReadStatusTracker.MATCHED_FIRST | ReadStatusTracker.MATCHED_SECOND;
  private static final int NUMBER_OF_THREAD_LOCKS = 1 << 16;
  private static final int THREAD_LOCK_MASK = NUMBER_OF_THREAD_LOCKS - 1;

  /* Array of locks for multiple threads */
  protected final Object[] mThreadLocks;

  protected final boolean mOutputUnmapped;
  protected final ReadBlocker mFreqBlockerLeft;
  private final ArrayList<Pair<HashingRegion, File>> mChildren;
  private final boolean mOutputSam;
  private int mChildCount = 0;

  /**
   * Create a new {@link UnfilteredSingleEndOutputProcessor}
   * @param param {@link NgsParams}
   * @param stats map to put statistics into
   * @param outputUnmapped  true if unmapped should be output
   * @throws IOException If an IO error
   */
  public UnfilteredSingleEndOutputProcessor(NgsParams param, MapStatistics stats, boolean outputUnmapped) throws IOException {
    super(param, stats, false, param.outputParams().sam() || param.outputParams().bam());
    final int numSequences = (int) param.buildFirstParams().numberSequences();
    mOutputUnmapped = outputUnmapped;
    mOutputSam = param.outputParams().sam() || mParams.outputParams().bam();
    mFreqBlockerLeft = new ReadBlocker(numSequences, param.readFreqThreshold(), "left hits");
    mChildren = new ArrayList<>();
    mThreadLocks = new Object[NUMBER_OF_THREAD_LOCKS];
    for (int i = 0; i < mThreadLocks.length; ++i) {
      mThreadLocks[i] = new Object();
    }
  }

  @Override
  public void finish() throws IOException {
    sortRegions();
    if (mSingleThreadChild != null) {
      synchronized (mSingleThreadChild) {
        mSingleThreadChild.threadFinish();
      }
    }
    // Merge the child output files into one
    mChildren.sort(new TopNPairedEndOutputProcessorSync.RegionFileComparator());

    final File[] outputFiles = new File[mChildren.size()];
    for (int i = 0; i < outputFiles.length; ++i) {
      outputFiles[i] = mChildren.get(i).getB();
    }
    final FilterConcatIntermediateFiles alignmentsIntFiles;
    if (mOutputSam) {
      final long readIdOffset = Math.max(mParams.buildFirstParams().readerRestriction().getStart(), 0);
      final File outFile;
      if (mParams.outputParams().bam()) {
        outFile = mParams.outputParams().resultStreamHandler().file(NgsOutputParams.ALIGNMENTS_BAM_FILE_NAME);
      } else {
        outFile = mParams.outputParams().isCompressOutput()
            ? mParams.outputParams().resultStreamHandler().file(NgsOutputParams.ALIGNMENTS_SAM_FILE_NAME + FileUtils.GZ_SUFFIX)
                : mParams.outputParams().resultStreamHandler().file(NgsOutputParams.ALIGNMENTS_SAM_FILE_NAME);
      }
      alignmentsIntFiles = new SansFilterConcat(mParams, mUnmappedTracker, readIdOffset).filterConcat(outputFiles, outFile, mSharedResources.getHeader(), mParams.outputParams());
      //SamSingleEndOutputProcessor.indexSamFile(mParams, outFile);
    } else {
      alignmentsIntFiles = null;
      for (final File f : outputFiles) {
        if (!f.delete()) {
          throw new IOException("Could not delete temporary file: " + f);
        }
      }
    }
    final FilterConcatIntermediateFiles unmappedIntFiles;
    if (mOutputUnmapped) {
      unmappedIntFiles = writeUnmapped();
    } else {
      unmappedIntFiles = null;
    }
    if (mParams.outputParams().unify() && mOutputSam) {
      whizBangUnify(unmappedIntFiles, alignmentsIntFiles);
    }
    super.finish();
    mUnmappedTracker.calculateStatistics(false, true);
  }

  private FilterConcatIntermediateFiles writeUnmapped() throws IOException {
    return writeUnmapped(!mParams.outputParams().unify(), !mOutputSam, true);
  }

  @Override
  protected FilterConcatIntermediateFiles filterConcatNonMated(MapQScoringReadBlocker blockerLeft, MapQScoringReadBlocker blockerRight, File[] tempFiles, SingleEndTopRandomImplementation.HitRecord[] hitsToKeep, NamesInterface templateNames, File outFile) {
    throw new UnsupportedOperationException();
  }
  @Override
  public synchronized OutputProcessor threadClone(HashingRegion region) throws IOException {
    super.threadClone(region);
    final int currentChild = mChildCount++;
    // Check that output directory exists
    final File dir = mParams.outputParams().directory();
    if (!dir.exists() && !dir.mkdirs()) {
      throw new IOException("unable to create directory: " + dir);
    }
    final File out = TopNPairedEndOutputProcessorSync.determineTempFile(mParams, currentChild);
    final OutputStream outStream;
    if (mOutputSam) {
      outStream = FileUtils.createOutputStream(out, true);
    } else {
      outStream = NullStreamUtils.getNullOutputStream();
    }
    final SingleEndTempFileWriter samse = new SingleEndTempFileWriter(mParams,  mUnmappedTracker, mSharedResources);
    samse.initialiseAlignments(outStream, null);
    samse.setClipRegion(region);
    final com.rtg.ngs.UnfilteredSingleEndOutputProcessor.InnerUnfilteredSingleEndOutputProcessor inner = new com.rtg.ngs.UnfilteredSingleEndOutputProcessor.InnerUnfilteredSingleEndOutputProcessor(this, samse);
    mChildren.add(new Pair<>(region, out));
    Diagnostic.developerLog("InnerUnfilteredOutputProcessor created for region:" + region);
    return inner;
  }

  @Override
  public void threadFinish() throws IOException {
    try {
      finish();
    } finally {
      close();
    }
  }
  @Override
  public void close() throws IOException {
    mFreqBlockerLeft.close();
    mSharedResources.close();
  }

  @Override
  public String toString() {
    return "SamUnfilteredOutputProcessor" + StringUtils.LS;
  }

  @Override
  public synchronized void process(final long templateId, final String frame, final int readId, final int tStart, final int score, final int scoreIndel) throws IOException {
    if (mChildCount == 0) {
      mSingleThreadChild = (com.rtg.ngs.UnfilteredSingleEndOutputProcessor.InnerUnfilteredSingleEndOutputProcessor) threadClone(HashingRegion.NONE);
    }
    mSingleThreadChild.process(templateId, frame, readId, tStart, score, scoreIndel);
  }

  com.rtg.ngs.UnfilteredSingleEndOutputProcessor.InnerUnfilteredSingleEndOutputProcessor mSingleThreadChild = null;

  protected static final class InnerUnfilteredSingleEndOutputProcessor implements OutputProcessor {
    private long mTemplateId = -1;
    private final SingleEndTempFileWriter mTempWriter;
    private final Object[] mThreadLocks;
    private final ReadBlocker mFreqBlockerLeft;
    private final ReadStatusTracker mUnmappedTracker;

    public InnerUnfilteredSingleEndOutputProcessor(UnfilteredSingleEndOutputProcessor parent, SingleEndTempFileWriter tempWriter) {
      mTempWriter = tempWriter;
      mThreadLocks = parent.mThreadLocks;
      mFreqBlockerLeft = parent.mFreqBlockerLeft;
      mUnmappedTracker = parent.mUnmappedTracker;
    }

    @Override
    public void threadFinish() throws IOException {
      try {
        finish();
      } finally {
        close();
      }
    }

    /** Not supported in thread local version */
    @Override
    public OutputProcessor threadClone(final HashingRegion region) {
      throw new UnsupportedOperationException("Should never get called.");
    }

    @Override
    public void finish() {
      Diagnostic.developerLog("Child finish");
    }

    @Override
    public void process(final long templateId, final String frame, final int readId, final int tStart, final int score, final int scoreIndel) throws IOException {
      final boolean bFrame = "R".equals(frame);
      if (templateId != mTemplateId) {
        mTemplateId = templateId;
        mTempWriter.nextTemplateId(templateId);
      }
      final BinaryTempFileRecord record;
      final ReadBlocker blocker = mFreqBlockerLeft;
      synchronized (mThreadLocks[readId & THREAD_LOCK_MASK]) {
        if (blocker.isBlocked(readId)) {
          return;
        }
        record = mTempWriter.alignmentResultUnfiltered(readId, bFrame, tStart);
        if (record != null) {
          blocker.increment(readId);
        }
        mUnmappedTracker.addStatus(readId, MATCHED);
      }
      //if (record != null && ReelTwoLicense.isDeveloper()) {
        //If you want to see score indel flags in SAM files, uncomment these lines.
        //record.setAttribute("XI", scoreIndel);
      //}
    }

    @Override
    public void close() throws IOException {
      mTempWriter.close();
    }
  }
}

