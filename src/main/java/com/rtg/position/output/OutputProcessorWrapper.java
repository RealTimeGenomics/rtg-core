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
package com.rtg.position.output;

import java.io.IOException;

import com.reeltwo.jumble.annotations.TestClass;
import com.rtg.index.hash.ngs.OutputProcessor;
import com.rtg.launcher.HashingRegion;
import com.rtg.launcher.ISequenceParams;
import com.rtg.mode.Frame;
import com.rtg.reader.SequencesReader;

/**
 * Wrapper for long reads. Assumes the we build on the reads.
 */
@TestClass(value = {"com.rtg.ngs.NgsLongTest"})
public class OutputProcessorWrapper implements PositionWriter {

  private final OutputProcessor mEnclosed;
  private final int[] mTemplateLengths;
  private final int mStart;
  private final SequencesReader mReader;
  private final SequencesReader mReader2;
  private final boolean mPairedEnd;
  private final HashingRegion mRegion;

  /**
   * Constructor
   * @param proc processor to wrap
   * @param reads for reads
   * @param readsSecond for second reads
   * @param template for template
   * @throws IOException If an I/O error occurs
   */
  public OutputProcessorWrapper(final OutputProcessor proc, final ISequenceParams reads, final ISequenceParams readsSecond, final ISequenceParams template) throws IOException {
    mEnclosed = proc;
    //mEnclosed = new TopNOutputProcessor(params.ngsOutputParams(), params.ngsOutputParams().outStream(), new NullMarcher(), readsReader, templateReader);
    final SequencesReader readsReader = reads.reader();
    final SequencesReader templateReader = template.reader();
    if (template.region() != HashingRegion.NONE) {
      mTemplateLengths = templateReader.sequenceLengths(template.region().getStart(), template.region().getExclusiveEndId());
      mStart = (int) template.region().getStart();
    } else {
      mTemplateLengths = templateReader.sequenceLengths(0, template.reader().numberSequences());
      mStart = 0;
    }
    mRegion = template.region();
    mReader = readsReader.copy();
    mPairedEnd = readsSecond != null;
    if (mPairedEnd) {
      mReader2 = readsSecond.reader().copy();
    } else {
      mReader2 = null;
    }
  }

  /**
   * @return the enclosed output processor
   */
  public OutputProcessor getEnclosed() {
    return mEnclosed;
  }


  /**
   * Not important
   * @param queryId ignored
   */
  @Override
  public void endQuery(final int queryId) {
  }

  @Override
  public void endAll() {
    //apparently writeResults called in Ngs
  }

  /**
   * Not important
   * @return 0.0
   */
  @Override
  public double score() {
    return 0.0;
  }

  @Override
  public void write(final AbstractGappedRegion<?> region, final int templateSequenceId, final Frame queryFrame, final int queryLength, final int queryEffectiveLength) throws IOException {
    //query = template, subject = reads
    final SurrogateRegion surr = region.surrogateOutput(null, templateSequenceId, queryFrame, queryLength, queryEffectiveLength);
    if (!surr.scoreAllowed()) {
      return;
    }
    final int readId = region.mBuildSeqId / region.mSubjectFrames.length;
    //TODO find out mismatch score
    final double score = surr.score();
    final long templateLength = mTemplateLengths[templateSequenceId - mStart];
    final int templatePos;
    final int templateEnd;

    final int queryHitLength = region.mQueryEnd - region.mQueryStart;
    final int buildHitLength = region.mBuildEnd - region.mBuildStart;
    final int indelDiff = queryHitLength - buildHitLength;
    final String frameOut;
    final int lengthFromReader;
    if (mPairedEnd) {
      if ((readId & 1) == 0) {
        lengthFromReader = mReader.length(readId / 2);
      } else {
        lengthFromReader = mReader2.length(readId / 2);
      }
    } else {
      lengthFromReader = mReader.length(readId);
    }
    final int readLength;
    if (region.mSubjectFrames.length > 1) {
      final int readFrameIndex = region.mBuildSeqId % region.mSubjectFrames.length;
      final Frame readFrame = region.mSubjectFrames[readFrameIndex];
      frameOut = readFrame.display();
      readLength = lengthFromReader / 3;
    } else {
      frameOut = queryFrame.isForward() ? "F" : "R";
      readLength = lengthFromReader;
    }

    //if (queryFrame.isForward()) {
      templatePos = region.mQueryStart - region.mBuildStart;
      templateEnd = templatePos + readLength + indelDiff;
    //} else {
    //  templatePos = (int) templateLength - (region.mQueryStart + readLength + indelDiff - region.mBuildStart);
    //  templateEnd = templatePos + readLength + indelDiff;
    //}
    final long clipStart;
    if (templateSequenceId == mRegion.getStart() && mRegion.getStartClipPosition() != HashingRegion.MISSING && mRegion.getStartPaddedPosition() != 0) {
      clipStart = mRegion.getStartPaddedPosition();
    } else {
      clipStart = Integer.MIN_VALUE;
    }
    final long clipEnd;
    if (templateSequenceId == mRegion.getEnd() && mRegion.getEndClipPosition() != HashingRegion.MISSING && mRegion.getEndPaddedPosition() != templateLength) {
      clipEnd = mRegion.getEndPaddedPosition();
    } else {
      clipEnd = Integer.MAX_VALUE;
    }

    /*
    final StringBuilder sb = new StringBuilder();
    sb.append(String.format("XXXX %d %s %d %d %.2f", templateSequenceId, frameOut, readId, templatePos, score));
    sb.append(" (R: " + region + ")");
    sb.append(String.format(" [S: %d %.2f %s]", surr.queryId(), surr.score(), surr.scoreAllowed() ? "true" : "false"));
    */
    if (templatePos >= clipStart && templateEnd <= clipEnd) {
      final int newScore = (int) score; //clipping? max 63? huh?
      //sb.append(" notclipped");
      mEnclosed.process(templateSequenceId, frameOut, readId, templatePos, 0, newScore);
    }
    //System.err.println(sb.toString());
  }

  /**
   * Return a clone of this suitable for use in the reverse frame
   * @return the clone
   */
  @Override
  public PositionWriter reverseClone() {
    return this;
  }


}
