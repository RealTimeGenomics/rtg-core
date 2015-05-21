/*
 * Copyright (c) 2015. Real Time Genomics Limited.
 *
 * Use of this source code is bound by the Real Time Genomics Limited Software Licence Agreement
 * for Academic Non-commercial Research Purposes only.
 *
 * If you did not receive a license accompanying this file, a copy must first be obtained by email
 * from support@realtimegenomics.com.  On downloading, using and/or continuing to use this source
 * code you accept the terms of that license agreement and any amendments to those terms that may
 * be made from time to time by Real Time Genomics Limited.
 */
package com.rtg.reader;

import static com.rtg.launcher.CommonFlags.STDIO_NAME;

import java.io.File;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.util.Locale;

import org.apache.commons.lang.ArrayUtils;

import com.reeltwo.jumble.annotations.TestClass;
import com.rtg.launcher.CommonFlags;
import com.rtg.util.diagnostic.ErrorType;
import com.rtg.util.diagnostic.NoTalkbackSlimException;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.LineWriter;

/**
 * Wrapper for writing FASTA that can handle single end or paired end data
 */
@TestClass("com.rtg.reader.Sdf2FastaTest")
public class FastaWriterWrapper implements WriterWrapper {

  private static final String[] EXTS = {".fasta", ".fa"};

  private final boolean mIsPaired;
  private final SdfReaderWrapper mReader;
  protected final int mLineLength; // Maximum residues per line -- 0 denotes infinite line length
  protected final boolean mRename; // Replace sequence names with their SDF sequence ID
  protected final boolean mHasNames;
  protected final byte[] mCodeToBytes;

  private final LineWriter mLeft;
  private final LineWriter mRight;
  private final LineWriter mSingle;


  /**
   * Convenience wrapper for writing.
   * @param baseOutput base output file name.
   * @param reader the reader that this writer is writing from.
   * @param lineLength the maximum line length, 0 means no bound.
   * @param rename if true, rename sequences to their sequence id
   * @param gzip if true, compress the output.
   * @throws IOException if there is a problem constructing the writer.
   */
  public FastaWriterWrapper(File baseOutput, SdfReaderWrapper reader, int lineLength, boolean rename, boolean gzip) throws IOException {
    this(baseOutput, reader, lineLength, rename, gzip, EXTS);
  }

  protected FastaWriterWrapper(File baseOutput, SdfReaderWrapper reader, int lineLength, boolean rename, boolean gzip, String[] extensions) throws IOException {
    assert reader != null;
    assert extensions.length > 0;
    mReader = reader;
    mIsPaired = reader.isPaired();
    mHasNames = reader.hasNames();
    final long maxLength = reader.maxLength();
    if (maxLength > Integer.MAX_VALUE) {
      throw new NoTalkbackSlimException(ErrorType.SEQUENCE_LENGTH_ERROR);
    }
    mCodeToBytes = SdfSubseq.getByteMapping(reader.type(), false);
    mLineLength = lineLength;
    mRename = rename;

    // Strip any existing GZ or FASTA like suffix from the file name, but remember the primary type
    String output = baseOutput.toString();
    if (FileUtils.isGzipFilename(output)) {
      output = output.substring(0, output.length() - FileUtils.GZ_SUFFIX.length());
    }
    String ext = FileUtils.getExtension(output);
    if (ArrayUtils.contains(extensions, ext.toLowerCase(Locale.getDefault()))) {
      output = output.substring(0, output.lastIndexOf('.'));
    } else {
      ext = extensions[0];
    }

    if (mIsPaired) {
      if (CommonFlags.isStdio(output)) {
        throw new NoTalkbackSlimException("Sending paired-end data to stdout is not supported.");
      }
      mSingle = null;
      mLeft = getStream(output + "_1" + ext, gzip);
      mRight = getStream(output + "_2" + ext, gzip);
    } else {
      mLeft = null;
      mRight = null;
      mSingle = getStream(CommonFlags.isStdio(output) ? STDIO_NAME : (output + ext), gzip);
    }
  }

  static LineWriter getStream(final String name, boolean gzip) throws IOException {
    if (CommonFlags.isStdio(name)) {
      return new LineWriter(new OutputStreamWriter(FileUtils.getStdoutAsOutputStream()));
    }
    return new LineWriter(new OutputStreamWriter(FileUtils.createOutputStream(gzip ? new File(name + FileUtils.GZ_SUFFIX) : new File(name), gzip)));
  }

  @Override
  public void writeSequence(long seqId, byte[] dataBuffer, byte[] qualityBuffer) throws IllegalStateException, IOException {
    if (mIsPaired) {
      writeSequence(mReader.left(), seqId, mLeft, dataBuffer, qualityBuffer);
      writeSequence(mReader.right(), seqId, mRight, dataBuffer, qualityBuffer);
    } else {
      writeSequence(mReader.single(), seqId, mSingle, dataBuffer, qualityBuffer);
    }
  }

  protected void writeSequence(SequencesReader reader, long seqId, LineWriter writer, byte[] dataBuffer, byte[] qualityBuffer) throws IllegalArgumentException, IllegalStateException, IOException {
    final int length = reader.read(seqId, dataBuffer);
    for (int i = 0; i < length; i++) {
      dataBuffer[i] = mCodeToBytes[dataBuffer[i]];
    }
    final String name = !mRename && mHasNames ? reader.fullName(seqId) : ("" + seqId);
    writer.writeln(">" + name);
    if (mLineLength == 0) {
      writer.writeln(new String(dataBuffer, 0, length));
    } else {
      for (long k = 0; k < length; k += mLineLength) {
        writer.writeln(new String(dataBuffer, (int) k, Math.min(mLineLength, length - (int) k)));
      }
    }
  }

  @Override
  public void close() throws IOException {
    if (mLeft != null) {
      mLeft.close();
    }
    if (mRight != null) {
      mRight.close();
    }
    if (mSingle != null) {
      mSingle.close();
    }
  }
}

