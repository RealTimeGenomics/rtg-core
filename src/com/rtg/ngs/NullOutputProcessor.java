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

import java.util.ArrayList;

import com.rtg.index.hash.ngs.OutputProcessor;
import com.rtg.launcher.HashingRegion;
import com.rtg.util.diagnostic.Diagnostic;

/**
 * Does absolutely nothing but counting of hits, to allow profiling of index searching.
 *
 */
public class NullOutputProcessor implements OutputProcessor {

  private long mHitCount = 0;
  private final ArrayList<NullOutputProcessor> mChildren = new ArrayList<>();

  /**
   * Does nothing
   */
  public NullOutputProcessor() {
  }

  @Override
  public void process(final long templateId, final String frame, final int readId, final int tStart, final int score, final int scoreIndel) {
    mHitCount++;
    /*
    System.out.println(""
                + templateId + "\t"
                + frame + "\t"
                + readId + "\t"
                + tStart + "\t"
                + score + "\t"
                + scoreIndel
                );
     */
  }

  @Override
  public void finish() {
    long totalHits = mHitCount; // Single thread case?
    for (NullOutputProcessor o : mChildren) {
      totalHits += o.mHitCount;
    }
    Diagnostic.developerLog("NullOutputProcessor total number of hits:" + totalHits);
  }

  /**
   */
  @Override
  public OutputProcessor threadClone(HashingRegion region) {
    final NullOutputProcessor inner = new NullOutputProcessor();
    mChildren.add(inner);
    Diagnostic.developerLog("NullOutputProcessor created for region:" + region);
    return inner;
  }

  @Override
  public void threadFinish() {
  }

  /**
   */
  @Override
  public void close() {
  }

  @Override
  public String toString() {
    return "NullOutputProcessor";
  }
}

