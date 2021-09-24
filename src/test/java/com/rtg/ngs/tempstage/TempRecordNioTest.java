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
package com.rtg.ngs.tempstage;

import java.io.InputStream;
import java.io.OutputStream;

/**
 */
public class TempRecordNioTest extends AbstractTempRecordReadWriteTest {

  @Override
  public TempRecordWriter getWriter(OutputStream out, TempRecordReader.RecordFactory fact) {
    return new TempRecordWriterNio(out);
  }
  @Override
  public TempRecordReader getReader(InputStream in, TempRecordReader.RecordFactory fact) {
    return new TempRecordReaderNio(in, fact);
  }
}
