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
package com.rtg.variant.cnv;

import static com.rtg.util.StringUtils.TAB;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;

import com.rtg.mode.DNA;
import com.rtg.mode.SequenceType;
import com.rtg.reader.SequencesReader;
import com.rtg.reader.SequencesReaderFactory;
import com.rtg.util.StringUtils;

/**
 * Detects regions of N's in the genome and reports them. Currently only supports DNA.
 * <p>
 * Output is in form:<br>
 * <code>seqId oneBaseStartPos inclusiveOneBasedEndPos</code>
 */
public final class NBlockDetector {

  private static final int MAX_UNKNOWN = DNA.N.ordinal();
  private static final String LS = StringUtils.LS;

  private NBlockDetector() {
  }

  /**
   * Detect all N blocks in provided sequence reader larger than the <code>blockSize</code>
   * @param sequence sequence data
   * @param blockSize minimum block size for output
   * @param output destination for output
   * @throws IOException IO Exceptions happen sometimes
   */
  public static void detectNs(SequencesReader sequence, int blockSize, OutputStream output) throws IOException {
    if (sequence.type() == SequenceType.DNA) {
      final byte[] buf = new byte[(int) sequence.maxLength()];
      for (long i = 0; i < sequence.numberSequences(); ++i) {
        final int length = sequence.read(i, buf);
        detectNs(sequence.name(i), buf, 0, length, blockSize, output);
      }
    } else {
      throw new UnsupportedOperationException("Only support DNA sequences");
    }
  }

  /**
   * Detect all N blocks in provided sequence reader larger than the <code>blockSize</code>
   * @param name Sequence name
   * @param seq sequence data
   * @param start start position in sequence data. zero based
   * @param end end position in sequence data (exclusive). zero based
   * @param blockSize minimum block size for output
   * @param output destination for output
   * @throws IOException IO Exceptions happen sometimes
   */
  public static void detectNs(String name, byte[] seq, int start, int end, int blockSize, OutputStream output) throws IOException {
    int nStart = 0;
    boolean nRegion = false;
    for (int i = start; i < end; ++i) {
      final boolean n = seq[i] <= MAX_UNKNOWN;
      if (!nRegion && n) {
        nStart = i;
        nRegion = true;
      } else if (nRegion && !n) {
        if (i - nStart >= blockSize) {
          output.write((name + TAB + (nStart + 1) + TAB + i + LS).getBytes());
        }
        nRegion = false;
      }
    }
    if (nRegion) {
      if (end - nStart >= blockSize) {
        output.write((name + TAB + (nStart + 1) + TAB + end + LS).getBytes());
      }
    }
  }

  /**
   * Command line entry point
   * @param args arguments
   * @throws IOException If an IO exception happens
   */
  public static void main(String[] args) throws IOException {
    if (args.length < 3) {
      System.err.println("Usage: NBlockDetector sdf-dir output-file blockLength");
      return;
    }
    final File dir = new File(args[0]);
    final File output = new File(args[1]);
    final int blockSize = Integer.parseInt(args[2]);
    try (FileOutputStream outputStream = new FileOutputStream(output);
      SequencesReader reader = SequencesReaderFactory.createDefaultSequencesReader(dir)) {
      detectNs(reader, blockSize, outputStream);
    }
  }
}
