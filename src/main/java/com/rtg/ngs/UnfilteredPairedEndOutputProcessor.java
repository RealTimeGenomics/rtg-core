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

import static com.rtg.index.hash.ngs.ReadDecoder.PAIRED_END;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;

import com.rtg.index.hash.ngs.OutputProcessor;
import com.rtg.launcher.HashingRegion;
import com.rtg.ngs.blocking.MapQScoringReadBlocker;
import com.rtg.ngs.blocking.ReadBlocker;
import com.rtg.ngs.blocking.ReadBlockerSync;
import com.rtg.ngs.tempstage.UnfilteredTempFileWriter;
import com.rtg.pairedend.ReadStatusListener;
import com.rtg.pairedend.UnfilteredSlidingWindowCollector;
import com.rtg.reader.NamesInterface;
import com.rtg.util.NullStreamUtils;
import com.rtg.util.Pair;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.io.FileUtils;

/**
 * Thread safe version of super class
 */
public class UnfilteredPairedEndOutputProcessor extends AbstractSdfOutputProcessor {

  //make dynamic
  final List<Pair<HashingRegion, File>> mChildren;

  private int mChildCount = 0;

  final boolean mOutputUnmapped;
  private final boolean mOutputSam;

  private final ReadBlocker mFreqBlockerLeft;
  private final ReadBlocker mFreqBlockerRight;

  /**
   * Construct an output processor
   * @param param parameters
   * @param stats map to put statistics into
   * @param outputUnmapped  true if unmapped should be output
   * @throws IOException if sam has problems setting up.
   */
  public UnfilteredPairedEndOutputProcessor(NgsParams param, MapStatistics stats, boolean outputUnmapped) throws IOException {
    super(param, stats, param.paired(), param.outputParams().sam() || param.outputParams().bam());
    final int sequences = (int) param.buildFirstParams().numberSequences();
    mChildren = new ArrayList<>();
    mOutputUnmapped = outputUnmapped;
    mOutputSam = param.outputParams().sam() || param.outputParams().bam();

    // These blockers are for counting hits per read per side, for high frequency filtering purposes
    mFreqBlockerLeft = new ReadBlockerSync(sequences, param.readFreqThreshold(), "left hits");
    mFreqBlockerRight = new ReadBlockerSync(sequences, param.readFreqThreshold(), "right hits");
  }

  private static synchronized void createDir(File dir) throws IOException {
    if (!dir.exists() && !dir.mkdirs()) {
      throw new IOException("unable to create directory: " + dir);
    }
  }

  @Override
  public synchronized OutputProcessor threadClone(HashingRegion region) throws IOException {
    super.threadClone(region);
    final int currentChild = mChildCount++;

    // Check that output directory exists
    final File dir = mParams.outputParams().directory();
    createDir(dir);
    final File out = TopNPairedEndOutputProcessorSync.determineTempFile(mParams, currentChild);
    final OutputStream outStream;
    if (mOutputSam) {
      outStream = FileUtils.createOutputStream(out, true);
    } else {
      outStream = NullStreamUtils.getNullOutputStream();
    }
    final UnfilteredTempFileWriter sam = new UnfilteredTempFileWriter(mUnmappedTracker, mSharedResources, mParams, mFreqBlockerLeft, mFreqBlockerRight);
    sam.initialiseUnmated(outStream);
    sam.setClipRegion(region);
    final UnfilteredSlidingWindowCollector collector = new UnfilteredSlidingWindowCollector(mParams.maxFragmentLength(), mParams.minFragmentLength(), mParams.pairOrientation(), sam,
        mSharedResources, mParams.outputParams().calibrateRegions());
    final InnerPairedEndOutputProcessor ipeop = new InnerPairedEndOutputProcessor(sam, collector, mFreqBlockerLeft, mFreqBlockerRight, mUnmappedTracker, region, out);
    mChildren.add(new Pair<>(ipeop.getRegion(), ipeop.getFile()));
    Diagnostic.developerLog("InnerPairedEndOutputProcessor:" + ipeop + " region:" + region);
    return ipeop;
  }

  @Override
  protected FilterConcatIntermediateFiles filterConcatNonMated(MapQScoringReadBlocker blockerLeft, MapQScoringReadBlocker blockerRight, File[] tempFiles, SingleEndTopRandomImplementation.HitRecord[] hitsToKeep, NamesInterface templateNames, File outFile) {
    throw new UnsupportedOperationException("Unpossible");
  }

  @Override
  public void threadFinish() throws IOException {
    try {
      finish();
    } finally {
      close();
    }
  }

  static class RegionFileComparator implements Comparator<Pair<HashingRegion, File>>, Serializable {
    @Override
    public int compare(Pair<HashingRegion, File> first, Pair<HashingRegion, File> second) {
      return first.getA().compareTo(second.getA());
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
    mChildren.sort(new RegionFileComparator());

    final File[] outputFiles = new File[mChildren.size()];
    for (int i = 0; i < outputFiles.length; ++i) {
      outputFiles[i] = mChildren.get(i).getB();
    }

    //merge results
    final FilterConcatIntermediateFiles alignmentsIntFiles;
    if (mOutputSam) {
      alignmentsIntFiles = writeAlignments(outputFiles);
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
    mUnmappedTracker.calculateStatistics(true, true);
  }

  private FilterConcatIntermediateFiles writeAlignments(File[] outputFiles) throws IOException {
    Diagnostic.userLog("Merging unmated results");
    final long readIdOffset = Math.max(mParams.buildFirstParams().readerRestriction().getStart(), 0);
    final File outFile;
    if (mParams.outputParams().bam()) {
      outFile = mParams.outputParams().resultStreamHandler().file(NgsOutputParams.ALIGNMENTS_BAM_FILE_NAME);
    } else {
      outFile = mParams.outputParams().isCompressOutput()
        ? mParams.outputParams().resultStreamHandler().file(NgsOutputParams.ALIGNMENTS_SAM_FILE_NAME + FileUtils.GZ_SUFFIX)
        : mParams.outputParams().resultStreamHandler().file(NgsOutputParams.ALIGNMENTS_SAM_FILE_NAME);
    }
    return new SansFilterConcat(mParams, mUnmappedTracker, readIdOffset).filterConcat(outputFiles, outFile, mSharedResources.getHeader(), mParams.outputParams());
  }

  private FilterConcatIntermediateFiles writeUnmapped() throws IOException {
    return writeUnmapped(!mParams.outputParams().unify(), !mOutputSam, true);
  }

  @Override
  public void close() throws IOException {
    mFreqBlockerLeft.close();
    mFreqBlockerRight.close();
    mSharedResources.close();
  }

  InnerPairedEndOutputProcessor mSingleThreadChild = null;

  @Override
  public synchronized void process(final long templateId, final String frame, final int readId, final int tStart, final int score, final int scoreIndel) throws IOException {
    if (mChildCount == 0) {
      mSingleThreadChild = (InnerPairedEndOutputProcessor) threadClone(HashingRegion.NONE);
    }
    mSingleThreadChild.process(templateId, frame, readId, tStart, score, scoreIndel);
    //throw new UnsupportedOperationException("Calls should be made on thread local clones");
  }

  static class InnerPairedEndOutputProcessor implements OutputProcessor {
    private final ReadStatusListener mListener;
    private final HashingRegion mRegion;
    private final File mFile;
    private final PairedEndOutputProcessor mPairedEndOutputProcessor;

    protected final ReadBlocker mFreqBlockerLeft;
    protected final ReadBlocker mFreqBlockerRight;
    private int mProcessCalls = 0;

    InnerPairedEndOutputProcessor(UnfilteredTempFileWriter sam, UnfilteredSlidingWindowCollector collector, ReadBlocker left, ReadBlocker right, ReadStatusListener listener, HashingRegion region, File file) {
      mPairedEndOutputProcessor = new PairedEndOutputProcessor(sam, collector);
      mListener = listener;
      mRegion = region;
      mFile = file;
      mFreqBlockerLeft = left;
      mFreqBlockerRight = right;
    }

    /**
     * @return the file this output processor is writing to.
     */
    public File getFile() {
      return mFile;
    }

    /**
     * @return the HashingRegion this output processor is processing
     */
    public HashingRegion getRegion() {
      return mRegion;
    }

    /** Not supported in thread local version */
    @Override
    public OutputProcessor threadClone(final HashingRegion region) {
      throw new UnsupportedOperationException("Not supported.");
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
      mPairedEndOutputProcessor.close();
    }

    @Override
    public void finish() throws IOException {
      mPairedEndOutputProcessor.finish();
      Diagnostic.developerLog("Child finish TopNPEOPChild process calls: " + mProcessCalls);
    }

    @Override
    public void process(final long templateId, final String frame, final int readId, final int tStart, final int score, final int scoreIndel) throws IOException {
      final int dReadId = PAIRED_END.decode(readId);
      final boolean isFirst = PAIRED_END.isFirst(readId);
      if (mFreqBlockerLeft.isBlocked(dReadId) && mFreqBlockerRight.isBlocked(dReadId)) {
        return;
      }
      ++mProcessCalls;
      final boolean reverse = frame.startsWith("R");
      mListener.addStatus(dReadId, isFirst ? ReadStatusTracker.MATCHED_FIRST : ReadStatusTracker.MATCHED_SECOND);
      mPairedEndOutputProcessor.process(templateId, reverse, readId, tStart);
    }
  }
}
