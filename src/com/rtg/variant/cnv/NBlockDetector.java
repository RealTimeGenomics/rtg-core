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
      for (long i = 0; i < sequence.numberSequences(); i++) {
        sequence.seek(i);
        final int length = sequence.readCurrent(buf);
        detectNs(sequence.currentName(), buf, 0, length, blockSize, output);
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
    for (int i = start; i < end; i++) {
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
