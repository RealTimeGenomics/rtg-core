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
