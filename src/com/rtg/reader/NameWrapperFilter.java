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

import java.io.IOException;
import java.util.Map;

import com.reeltwo.jumble.annotations.TestClass;
import com.rtg.util.diagnostic.Diagnostic;

/**
 * Facilitates transferring sequences by (short) name from a SdfReaderWrapper to a WriterWrapper.
 */
@TestClass("com.rtg.reader.SdfSubsetTest")
class NameWrapperFilter extends WrapperFilter {

  private final Map<String, Long> mNames;
  private final AbstractSdfWriter.SequenceNameHandler mHandler = new AbstractSdfWriter.SequenceNameHandler();

  NameWrapperFilter(SdfReaderWrapper reader, WriterWrapper writer) throws IOException {
    super(reader, writer);
    mNames = ReaderUtils.getSequenceNameMap(reader.isPaired() ? reader.left() : reader.single());
  }

  @Override
  protected void warnInvalidSequence(String seqid) {
    if (mWarnCount < 5) {
      Diagnostic.warning("Invalid sequence name " + seqid);
      mWarnCount++;
      if (mWarnCount == 5) {
        Diagnostic.warning("(Only the first 5 messages shown.)");
      }
    } else {
      Diagnostic.userLog("Invalid sequence name " + seqid);
    }
  }

  /**
   * Transfer an interpreted sequence or set of sequences from the reader to the writer.
   * This implementation interprets the specifier as a short sequence name.
   * @param seqRange the sequence name
   * @throws IOException if there was a problem during writing
   */
  @Override
  protected void transfer(String seqRange) throws IOException {
    final Long id = mNames.get(mHandler.handleSequenceName(seqRange).label());
    if (id != null) {
      transfer(id);
    } else {
      warnInvalidSequence(seqRange);
    }
  }
}
