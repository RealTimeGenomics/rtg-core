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
package com.rtg.simulation.reads;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.Writer;

import com.rtg.mode.DnaUtils;
import com.rtg.reader.FastaUtils;
import com.rtg.reader.SdfId;
import com.rtg.util.io.LineWriter;

/**
 * FASTQ style read simulator output
 */
public class FastqReadWriter implements ReadWriter {

  private Appendable mAppend;
  private final LineWriter mAppendLeft;
  private final LineWriter mAppendRight;
  private int mTotal = 0;
  private boolean mExpectLeft = true;

  /**
   * Constructor
   * @param append destination for output
   */
  public FastqReadWriter(Appendable append) {
    mAppend = append;
    mAppendLeft = null;
    mAppendRight = null;
  }

  /**
   * Constructor for use with two appendables (e.g. paired end <code>fastq</code> output)
   * @param fastqBaseFileName base name of a paired end <code>fastq</code> file to output to
   * @throws IOException if an exception occurs when instantiating writers
   */
  public FastqReadWriter(File fastqBaseFileName) throws IOException {
    final String fpath = fastqBaseFileName.getPath();
    final String base = fpath.substring(0, fpath.length() - 3);
    mAppendLeft = new LineWriter(new FileWriter(base + "_1.fq"));
    mAppendRight = new LineWriter(new FileWriter(base + "_2.fq"));
  }

  @Override
  public void identifyTemplateSet(SdfId... templateIds) {
    // Ignored
  }

  @Override
  public void identifyOriginalReference(SdfId referenceId) {
    // Ignored
  }

  @Override
  public void writeLeftRead(String name, byte[] data, byte[] qual, int length) throws IOException {
    if (!mExpectLeft) {
      throw new IllegalStateException();
    }
    if (mAppendLeft != null) {
      mAppend = mAppendLeft;
    }
    writeSequence(mTotal + " " + name + " 1", data, qual, length);
    mExpectLeft = !mExpectLeft;
  }

  @Override
  public void writeRightRead(String name, byte[] data, byte[] qual, int length) throws IOException {
    if (mExpectLeft) {
      throw new IllegalStateException();
    }
    if (mAppendRight != null) {
      mAppend = mAppendRight;
    }
    writeSequence(mTotal + " " + name + " 2", data, qual, length);
    mExpectLeft = !mExpectLeft;
    mTotal++;
  }

  @Override
  public void writeRead(String name, byte[] data, byte[] qual, int length) throws IOException {
    writeSequence(mTotal + " " + name, data, qual, length);
    mTotal++;
  }

  private void writeSequence(String name, byte[] data, byte[] qual, int length) throws IOException {
    mAppend.append("@");
    mAppend.append(name);
    mAppend.append("\n");
    mAppend.append(DnaUtils.bytesToSequenceIncCG(data, 0, length));
    mAppend.append("\n");
    mAppend.append("+");
    //mAppend.append(name);
    mAppend.append("\n");
    mAppend.append(FastaUtils.rawToAsciiString(qual, 0, length));
    mAppend.append("\n");
  }

  @Override
  @SuppressWarnings("try")
  public void close() throws IOException {
    try (Writer ignored = mAppendLeft; Writer ignored2 = mAppendRight) {
      // we want the sexy closing side effects
    }
  }


  @Override
  public int readsWritten() {
    return mTotal;
  }

}
