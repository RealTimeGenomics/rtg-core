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

import java.io.File;
import java.io.IOException;

import com.reeltwo.jumble.annotations.TestClass;
import com.rtg.mode.SequenceType;
import com.rtg.util.InvalidParamsException;
import com.rtg.util.StringUtils;
import com.rtg.util.diagnostic.ErrorType;
import com.rtg.util.io.LineWriter;

/**
 * Wrapper for writing FASTA that can handle single end or paired end data
 */
@TestClass("com.rtg.reader.Sdf2FastqTest")
public final class FastqWriterWrapper extends FastaWriterWrapper {

  private static final String[] EXTS = {".fastq", ".fq"};

  private String mDefaultQualities = null;

  /**
   * Convenience wrapper for writing.
   * @param baseOutput base output file name.
   * @param reader the reader that this writer is writing from.
   */
  public FastqWriterWrapper(File baseOutput, SdfReaderWrapper reader, int lineLength, boolean rename, boolean gzip, int def) throws IOException {
    super(baseOutput, reader, lineLength, rename, gzip, EXTS);
    if (reader.type() != SequenceType.DNA) {
      throw new InvalidParamsException(ErrorType.INFO_ERROR, "The input SDF contains protein data, which cannot be converted to FASTQ.");
    }
    if (!reader.hasQualityData()) {
      if (def >= (int) '!') {
        mDefaultQualities = StringUtils.getCharString((char) def, reader.maxLength());
      } else {
        throw new InvalidParamsException(ErrorType.INFO_ERROR, "The input SDF does not have quality data and no default was provided.");
      }
    }
  }

  @Override
  protected void writeSequence(SequencesReader reader, long seqId, LineWriter writer, byte[] dataBuffer, byte[] qualityBuffer) throws IllegalArgumentException, IllegalStateException, IOException {
    final int length = reader.read(seqId, dataBuffer);
    for (int i = 0; i < length; i++) {
      dataBuffer[i] = mCodeToBytes[dataBuffer[i]];
    }
    final String name = !mRename && mHasNames ? reader.fullName(seqId) : ("" + seqId);
    writer.writeln("@" + name);
    if (mLineLength == 0) {
      writer.writeln(new String(dataBuffer, 0, length));
      writer.writeln("+" + name);
      writer.writeln(getScore(reader, seqId, qualityBuffer));
    } else {
      for (long k = 0; k < length; k += mLineLength) {
        writer.writeln(new String(dataBuffer, (int) k, Math.min(mLineLength, length - (int) k)));
      }
      writer.writeln("+" + name);
      final String qual = getScore(reader, seqId, qualityBuffer);
      final int qualLen = qual.length();
      for (long i = 0; i < qualLen; i += mLineLength) {
        writer.writeln(qual.substring((int) i, (int) Math.min(i + mLineLength, qualLen)));
      }
    }
  }

  private String getScore(final SequencesReader read, long seqId, byte[] qualityBuffer) throws IOException {
    if (read.hasQualityData()) {
      final int length = read.readQuality(seqId, qualityBuffer);
      return FastaUtils.rawToAsciiString(qualityBuffer, length);
    } else {
      return mDefaultQualities.substring(0, read.length(seqId));
    }
  }

}

