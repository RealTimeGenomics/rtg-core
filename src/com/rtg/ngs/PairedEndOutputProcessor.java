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
