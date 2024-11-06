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

import java.io.IOException;
import java.util.Map.Entry;

import com.rtg.ngs.tempstage.AbstractTempFileWriter;
import com.rtg.pairedend.AbstractSlidingWindowCollector;
import com.rtg.util.diagnostic.Diagnostic;

/**
 * Output processor that understands simultaneous paired-end hit collection.
 *
 */
public class PairedEndOutputProcessor {

  protected AbstractTempFileWriter mSamWriter;
  private AbstractSlidingWindowCollector<?> mCollector;

  private long mCurrentTemplateId = Long.MAX_VALUE;

  /**
   * @param sam thing to close at end
   * @param collector where we send hits for further processing
   */
  public PairedEndOutputProcessor(final AbstractTempFileWriter sam, final AbstractSlidingWindowCollector<?> collector) {
    //TODO collectors generally are constructed with ASAWs, can we simplify this?
    mSamWriter = sam;
    mCollector = collector;
  }

  /**
   * Process the given results.
   * @param templateId id of template
   * @param reverse true if reverse complement hit
   * @param readId read number (0 based).
   * @param tStart template start (1 based).
   * @throws IOException if problems writing output
   */
  public void process(final long templateId, boolean reverse, final int readId, final int tStart) throws IOException {
    if (templateId != mCurrentTemplateId) {
      mCurrentTemplateId = templateId;
      mCollector.nextTemplateId(mCurrentTemplateId);
    }

    assert readId >= 0;
    final boolean first = PAIRED_END.isFirst(readId);
    final int lReadId = PAIRED_END.decode(readId);
    //System.out.println("templateId=" + templateId + " frame=" + frame + " readId=" + lReadId + " tStart=" + tStart + " score=" + score + " scoreIndel=" + scoreIndel + " first=" + first);
    assert lReadId < Integer.MAX_VALUE;
    mCollector.match(first, reverse, lReadId, tStart);
  }

  /**
   * Informs the output processor that hit processing is finished and it is time to do any
   * end of run processing/cleanup/summarizing
   * @throws IOException if problems when writing output.
   */
  public void finish() throws IOException {
    if (mCollector != null) {
      mCollector.nextTemplateId(Long.MAX_VALUE);
      Diagnostic.developerLog("Sliding window collector statistics");
      for (final Entry<Object, Object> e : mCollector.getStatistics().entrySet()) {
        Diagnostic.developerLog(e.getKey() + " = " + e.getValue());
      }
    }
    mCollector = null;
  }

  /**
   * Closes anything that needs to be closed.
   * @throws IOException if an I/O error occurs
   */
  public void close() throws IOException {
    if (mSamWriter != null) {
      mSamWriter.close();
    }
    mSamWriter = null;
  }
}
