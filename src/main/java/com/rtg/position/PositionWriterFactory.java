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
package com.rtg.position;

import java.io.IOException;

import com.rtg.index.hash.ngs.OutputProcessor;
import com.rtg.launcher.ISequenceParams;
import com.rtg.launcher.HashingRegion;
import com.rtg.position.output.OutputProcessorWrapper;
import com.rtg.position.output.PositionWriter;

/**
 * A factory to create PositionWriters.
 * Needed because <code>Ngs</code> needs to call <code>Position</code>
 * for long reads.
 */
public abstract class PositionWriterFactory {

  /** Used for long reads */
  public static class NgsPositionWriterFactory extends PositionWriterFactory {
    private final OutputProcessor mProcessor;
    /**
     * Constructor
     * @param proc output processor to wrap
     */
    public NgsPositionWriterFactory(final OutputProcessor proc) {
      mProcessor = proc;
    }

    @Override
    public PositionWriter makeNgs(final ISequenceParams build, ISequenceParams buildSecond, final ISequenceParams query, HashingRegion r) throws IOException {
      return new OutputProcessorWrapper(mProcessor.threadClone(r), build, buildSecond, query);
    }

    @Override
    public OutputProcessor makeProcessor(HashingRegion r) throws IOException {
      return mProcessor.threadClone(r);
    }
  }

  /**
   * @param build Sequences reader for build
   * @param buildSecond second reader for build, if paired
   * @param query sequences reader for query
   * @param r HashingRegion on which output processor will handle
   * @return a <code>PositionWriter</code> which handles all output.
   * @throws IOException If an I/O error occurs
   */
  public abstract PositionWriter makeNgs(ISequenceParams build, ISequenceParams buildSecond, ISequenceParams query, HashingRegion r) throws IOException;

  /**
   * Make an output processor for the current thread
   * @param r the region the current thread is processing
   * @return an output processor for the current thread
   * @throws IOException if an I/O error occurs
   */
  public abstract OutputProcessor makeProcessor(HashingRegion r) throws IOException;
}
