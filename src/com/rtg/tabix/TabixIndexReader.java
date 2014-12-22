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
package com.rtg.tabix;

import java.io.ByteArrayOutputStream;
import java.io.EOFException;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.Arrays;
import java.util.HashMap;

import com.rtg.util.gzip.GzipUtils;
import com.rtg.util.io.ByteArrayIOUtils;

/**
 * Reads a <code>Tabix</code> index to retrieve appropriate information
 */
public class TabixIndexReader extends AbstractIndexReader {

  //private static final int LINEAR_INDEX_SHIFT = 14;

  private final TabixIndexer.TabixOptions mOptions;

//  private final int COORDINATE_MASK = 0x10000;
//  private final int FORMAT_MASK = 0xFFFF;

  /** Size of fixed portion of header */
  private static final int FIXED_HEADER_SIZE = 36;
  private static final byte[] TBI_MAGIC = {'T', 'B', 'I', 1};

  /**
   * @param tabixFile <code>TABIX</code> file to open
   * @throws IOException if an IO error occurs
   */
  public TabixIndexReader(File tabixFile) throws IOException {
    super(tabixFile);
    if (!TabixIndexer.isBlockCompressed(tabixFile)) {
      throw new IOException("File: " + tabixFile.getPath() + " is not a valid TABIX index. (index file is not block compressed)");
    }
    final byte[] bufIn = new byte[10 * 1024];
    final ByteArrayOutputStream bufHelper = new ByteArrayOutputStream();
    try (InputStream is = GzipUtils.createGzipInputStream(new FileInputStream(tabixFile))) {
      int inLen;
      while ((inLen = is.read(bufIn, 0, bufIn.length)) != -1) {
        bufHelper.write(bufIn, 0, inLen);
      }
      final byte[] buf = bufHelper.toByteArray();
      final int len = buf.length;
      if (len < FIXED_HEADER_SIZE) {
        throw new EOFException("File: " + tabixFile.getPath() + " is not a valid TABIX index. (index does not have a complete header)");
      }
      for (int i = 0; i < TBI_MAGIC.length; i++) {
        if (TBI_MAGIC[i] != buf[i]) {
          throw new IOException("File: " + tabixFile.getPath() + " is not a valid TABIX index. (index does not have a valid header)");
        }
      }
      final int numberReferences = ByteArrayIOUtils.bytesToIntLittleEndian(buf, 4);
      final int format = ByteArrayIOUtils.bytesToIntLittleEndian(buf, 8);
//      mSam = (format & FORMAT_MASK) == 1;
//      mVcf = (format & FORMAT_MASK) == 2;
//      mGeneric = (format & FORMAT_MASK) == 0;
      final int seqCol = ByteArrayIOUtils.bytesToIntLittleEndian(buf, 12) - 1;
      final int begCol = ByteArrayIOUtils.bytesToIntLittleEndian(buf, 16) - 1;
      final int endCol = ByteArrayIOUtils.bytesToIntLittleEndian(buf, 20) - 1;
      final int meta = ByteArrayIOUtils.bytesToIntLittleEndian(buf, 24);
      final int skip = ByteArrayIOUtils.bytesToIntLittleEndian(buf, 28);
      mOptions = new TabixIndexer.TabixOptions(format, seqCol, begCol, endCol, meta, skip);
      final int nameLength = ByteArrayIOUtils.bytesToIntLittleEndian(buf, 32);
      if (buf.length < FIXED_HEADER_SIZE + nameLength) {
        throw new IOException("File: " + tabixFile.getPath() + " is not a valid TABIX index. (sequence names truncated)");
      }
      final byte[] nameBuf = Arrays.copyOfRange(buf, FIXED_HEADER_SIZE, FIXED_HEADER_SIZE + nameLength);
      mSequenceNames = new String[numberReferences];
      mSequenceLookup = new HashMap<>();
      int seqNo = 0;
      int start = 0;
      for (int i = 0; i < nameLength; i++) {
        if (nameBuf[i] == 0) {
          mSequenceNames[seqNo] = new String(nameBuf, start, i - start);
          mSequenceLookup.put(mSequenceNames[seqNo], seqNo);
          start = i + 1;
          seqNo++;
        }
      }
      seqNo = 0;
      mBinPositions = new long[numberReferences];
      mLinearIndexPositions = new long[numberReferences];
      long pos = FIXED_HEADER_SIZE + nameLength;
      for (; seqNo < numberReferences; seqNo++) {
        mBinPositions[seqNo] = pos;
        final int numBins = ByteArrayIOUtils.bytesToIntLittleEndian(buf, (int) pos);
        pos += 4;
        for (int i = 0; i < numBins; i++) {
          final int numChunks = ByteArrayIOUtils.bytesToIntLittleEndian(buf, (int) pos + 4);
          pos += 8;
          pos += numChunks * 16L;
        }
        mLinearIndexPositions[seqNo] = pos;
        final int numLinear = ByteArrayIOUtils.bytesToIntLittleEndian(buf, (int) pos);
        pos += 4;
        pos += numLinear * 8L;
      }
    }
  }

  /**
   * @return options defined in <code>TABIX</code> file
   */
  public TabixIndexer.TabixOptions getOptions() {
    return mOptions;
  }

  @Override
  public InputStream openIndexFile() throws IOException {
    return GzipUtils.createGzipInputStream(new FileInputStream(mIndexFile));
  }
}
