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

import com.rtg.launcher.HashingRegion;
import com.rtg.ngs.tempstage.AbstractTempFileWriter;
import com.rtg.util.IORunnable;
import com.rtg.util.ProgramState;
import com.rtg.util.StringUtils;
import com.rtg.util.diagnostic.Diagnostic;

/**
 * Used to handle the batching of multithreaded writing of unmated or
 * single-end alignments.
 *
 */
abstract class AbstractAlignmentWriterThread implements IORunnable {

  protected final AbstractTempFileWriter mSamWriter;
  private final MatchResult mResults;
  private final int mThreadNumber;
  private final AlignmentWorkload mWorkload;


  /**
   * Creates a new <code>AbstractAlignmentWriterThread</code> instance.
   *
   * @param writer an <code>AbstractSamAlignmentWriter</code> that will receive the results to write
   * @param results the set of all results to be aligned and written, which will be partitioned among the threads
   * @param threadNumber the number of this thread
   */
  AbstractAlignmentWriterThread(AbstractTempFileWriter writer, MatchResult results, long chunkStart, long chunkEnd, HashingRegion region, int threadNumber) {
    mSamWriter = writer;
    mResults = results;
    mThreadNumber = threadNumber;
    mWorkload = new AlignmentWorkload(results, chunkStart, chunkEnd, region, threadNumber);
  }

  @Override
  public void run() throws IOException {
    final HashingRegion region = mWorkload.toRegion();
    final String name = getName() + region;
    mSamWriter.setClipRegion(region);
    Diagnostic.userLog(name + " starting, hits:" + mWorkload.mChunkStart + "-" + mWorkload.mChunkEnd);
    try {
      int templateId = -1;
      for (long i = mWorkload.getChunkStart(); i > -1 && i < mResults.size() && i < mWorkload.getChunkEnd(); i++) {
        if ((i & 0xFFFL) == 0) {
          ProgramState.checkAbort();
        }
        if (templateId != mResults.getTemplateId(i)) {
          templateId = mResults.getTemplateId(i);
          mSamWriter.nextTemplateId(templateId);
        }
        handleResult(templateId, mResults.getEncodedReadId(i), mResults.getPosition(i), mResults.isReverse(i));
      }
    } finally {
      mSamWriter.close();
    }
    Diagnostic.userLog(name + " finished");
  }

  protected abstract void handleResult(int templateId, int encodedReadId, int position, boolean reverse) throws IOException;


  @Override
  public String toString() {
    return String.valueOf(mThreadNumber);
  }

  protected String getName() {
    return "Alignment Processing Thread ";
  }

  static final String LS = StringUtils.LS;

  public String makeString() {
    return mWorkload.makeString();
  }

  static class AlignmentWorkload {
    final MatchResult mResults;
    final int mThreadNumber;
    final int mThreadCount;
    final long mChunkStart;
    final long mChunkEnd;
    final HashingRegion mRegion;

    AlignmentWorkload(MatchResult results, long chunkStart, long chunkEnd, HashingRegion region, int threadNumber) {
      mResults = results;
      mChunkStart = chunkStart;
      mChunkEnd = chunkEnd;
      mRegion = region;
      mThreadNumber = threadNumber;
      mThreadCount = -1;
    }

    public AlignmentWorkload(int nt, int threadPadding, MatchResult results, int threadNumber) {
      mResults = results;
      mThreadNumber = threadNumber;
      mThreadCount = nt;
      if (results.size() == 0) {
        mChunkStart = -1;
        mChunkEnd = -1;

        mRegion = HashingRegion.NONE;
      } else {
        final long start = splitStart();
        final long end;
        if (threadNumber == nt - 1) {
          end = mResults.size();
        } else {
          end = splitStart() + splitSize();
        }

        if (start < mResults.size()) {
          final int startTemplate = mResults.getTemplateId(start);
          final int startPosition = start == 0 ? 0 : mResults.getPosition(start);
          final int endTemplate;
          final int endPosition;
          if (end < mResults.size()) {
            endTemplate = mResults.getTemplateId(end);
            endPosition = mResults.getPosition(end);
          } else {
            endTemplate = mResults.getTemplateId(mResults.size() - 1);
            endPosition = mResults.getPosition(mResults.size() - 1) + threadPadding;
          }

          mRegion = new HashingRegion(startTemplate, startPosition, endTemplate, endPosition, Math.max(startPosition - threadPadding, 0), endPosition + threadPadding);

          long chunkStart = start;
          while (chunkStart > 0 && mRegion.isInPaddedRange(mResults.getTemplateId(chunkStart - 1), mResults.getPosition(chunkStart - 1)) == 0) {
            chunkStart--;
          }
          long chunkEnd = end;
          while (chunkEnd < mResults.size() && mRegion.isInPaddedRange(mResults.getTemplateId(chunkEnd), mResults.getPosition(chunkEnd)) == 0) {
            chunkEnd++;
          }
          mChunkStart = chunkStart;
          mChunkEnd = chunkEnd;
        } else {
          // do something for stupid case
          mRegion = HashingRegion.NONE;
          mChunkStart = mResults.size();
          mChunkEnd = mResults.size();
        }
      }
    }

    public long getChunkStart() {
      return mChunkStart;
    }

    public long getChunkEnd() {
      return mChunkEnd;
    }

    public HashingRegion toRegion() {
      return mRegion;
    }
    private long splitSize() {
      return mResults.size() / mThreadCount > 0 ? mResults.size() / mThreadCount : 1L;
    }
    private long splitStart() {
      return splitSize() * mThreadNumber;
    }
    public String makeString() {
      //mSamWriter;
      return ("results.length= " + mResults.size() + LS) + "thread number= " + mThreadNumber + LS + "thread count= " + mThreadCount + LS + "chunk start= " + mChunkStart + LS + "chunk end= " + mChunkEnd + LS + "clip end position= " + mRegion.getEndClipPosition() + LS + "end id= " + mRegion.getEnd() + LS + "start id= " + mRegion.getStart() + LS + "clip start position= " + mRegion.getStartClipPosition() + LS;
    }
  }
}
