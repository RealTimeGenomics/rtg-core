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
   * @param chunkStart start position of this chunk of output
   * @param chunkEnd end position of this chunk of output
   * @param region contains region info regarding padding etc.
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
    try (final AbstractTempFileWriter out = mSamWriter) {
      out.setClipRegion(region);
      Diagnostic.userLog(name + " starting, hits:" + mWorkload.mChunkStart + "-" + mWorkload.mChunkEnd);
      int templateId = -1;
      for (long i = mWorkload.getChunkStart(); i > -1 && i < mResults.size() && i < mWorkload.getChunkEnd(); ++i) {
        if ((i & 0xFFFL) == 0) {
          ProgramState.checkAbort();
        }
        if (templateId != mResults.getTemplateId(i)) {
          templateId = mResults.getTemplateId(i);
          out.nextTemplateId(templateId);
        }
        handleResult(templateId, mResults.getEncodedReadId(i), mResults.getPosition(i), mResults.isReverse(i));
      }
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

    AlignmentWorkload(int nt, int threadPadding, MatchResult results, int threadNumber) {
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
            --chunkStart;
          }
          long chunkEnd = end;
          while (chunkEnd < mResults.size() && mRegion.isInPaddedRange(mResults.getTemplateId(chunkEnd), mResults.getPosition(chunkEnd)) == 0) {
            ++chunkEnd;
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
      return "results.length= " + mResults.size() + LS + "thread number= " + mThreadNumber + LS + "thread count= " + mThreadCount + LS + "chunk start= " + mChunkStart + LS + "chunk end= " + mChunkEnd + LS + "clip end position= " + mRegion.getEndClipPosition() + LS + "end id= " + mRegion.getEnd() + LS + "start id= " + mRegion.getStart() + LS + "clip start position= " + mRegion.getStartClipPosition() + LS;
    }
  }
}
