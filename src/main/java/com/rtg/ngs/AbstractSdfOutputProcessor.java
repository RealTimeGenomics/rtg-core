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

import com.rtg.reader.SdfReaderWrapper;
import com.rtg.reader.SdfWriterWrapper;
import com.rtg.util.diagnostic.Diagnostic;

/**
 * Abstract Output Processor which allows SDF writing
 */
public abstract class AbstractSdfOutputProcessor extends AbstractMapOutputProcessor {

  private static final String ALIGNMENTS_SDF_FILE = "alignments.sdf";
  /** Filename used for unmapped records */
  public static final String UNMAPPED_SDF_FILE = "unmapped.sdf";

  final NgsParams mParams;

  private final boolean mOutputSdf;


  /**
   * Create an abstract SDF output processor
   *
   * @param param ngs params
   * @param stats statistics container
   * @param paired whether this is a paired end read output processor
   * @param outputAlignments whether we are outputting SAM/BAM alignments
   * @throws IOException if an io exception occurs
   */
  public AbstractSdfOutputProcessor(NgsParams param, MapStatistics stats, boolean paired, boolean outputAlignments) throws IOException {
    super(param, stats, paired, outputAlignments);
    mParams = param;
    mOutputSdf = param.outputParams().sdf();
  }

  @Override
  public void finish() throws IOException {
    if (mOutputSdf) {
      Diagnostic.userLog("Starting SDF output");
      final SdfReaderWrapper reader = new SdfReaderWrapper(mSharedResources.firstReaderCopy(), mSharedResources.secondReaderCopy());
      final int maxLength = reader.maxLength();
      final byte[] dataBuffer = new byte[maxLength];
      final byte[] qualityBuffer = new byte[maxLength];
      try (SdfWriterWrapper alignments = new SdfWriterWrapper(mParams.file(ALIGNMENTS_SDF_FILE), reader, false);
           SdfWriterWrapper unmapped = new SdfWriterWrapper(mParams.file(UNMAPPED_SDF_FILE), reader, false)) {
        for (long seq = 0; seq < reader.numberSequences(); ++seq) {
          if (!mUnmappedTracker.getStatus((int) seq, ReadStatusTracker.UNMAPPED_FIRST)
            || !mUnmappedTracker.getStatus((int) seq, ReadStatusTracker.UNMAPPED_SECOND)) {
            alignments.writeSequence(seq, dataBuffer, qualityBuffer);
          } else {
            unmapped.writeSequence(seq, dataBuffer, qualityBuffer);
          }
        }
      }
      Diagnostic.userLog("Finished SDF output");
    }
  }
}
