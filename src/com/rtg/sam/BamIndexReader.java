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
package com.rtg.sam;

import java.io.EOFException;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.HashMap;

import com.rtg.tabix.AbstractIndexReader;
import com.rtg.util.io.ByteArrayIOUtils;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.IOUtils;

import net.sf.samtools.SAMSequenceDictionary;

/**
 * Used for reading file block compressed file coordinates from BAM indexes
 */
public class BamIndexReader extends AbstractIndexReader {

  /** Size of fixed portion of header */
  private static final int FIXED_HEADER_SIZE = 8;
  private static final byte[] BAI_MAGIC = {'B', 'A', 'I', 1};

  /**
   * @param indexFile BAM index file to open
   * @param dict sequence dictionary from BAM file
   * @throws IOException if an IO error occurs
   */
  public BamIndexReader(File indexFile, SAMSequenceDictionary dict) throws IOException {
    super(indexFile);
    try (InputStream is = new FileInputStream(indexFile)) {
      final byte[] buf = new byte[FIXED_HEADER_SIZE];
      final int len = IOUtils.readAmount(is, buf, 0, FIXED_HEADER_SIZE);
      if (len != FIXED_HEADER_SIZE) {
        throw new EOFException("File: " + indexFile.getPath() + " is not a valid BAM index. (index does not have a complete header)");
      }
      for (int i = 0; i < BAI_MAGIC.length; i++) {
        if (BAI_MAGIC[i] != buf[i]) {
          throw new IOException("File: " + indexFile.getPath() + " is not a valid BAI index. (index does not have a valid header)");
        }
      }
      final int numRefs = ByteArrayIOUtils.bytesToIntLittleEndian(buf, 4);
      if (numRefs != dict.size()) {
        throw new IllegalArgumentException("Index file: " + indexFile.getPath() + " does not contain the same number of sequences as given sequence dictionary");
      }

      mSequenceNames = new String[numRefs];
      mSequenceLookup = new HashMap<>();
      for (int i = 0; i < dict.size(); i++) {
        mSequenceNames[i] = dict.getSequence(i).getSequenceName();
        mSequenceLookup.put(mSequenceNames[i], i);
      }
      int seqNo = 0;
      mBinPositions = new long[numRefs];
      mLinearIndexPositions = new long[numRefs];
      long pos = FIXED_HEADER_SIZE;
      final byte[] bBuf = new byte[8];
      for (; seqNo < numRefs; seqNo++) {
        mBinPositions[seqNo] = pos;
        IOUtils.readFully(is, bBuf, 0, 4);
        pos += 4;
        final int numBins = ByteArrayIOUtils.bytesToIntLittleEndian(bBuf, 0);
        for (int i = 0; i < numBins; i++) {
          IOUtils.readFully(is, bBuf, 0, 8);
          pos += 8;
          final int numChunks = ByteArrayIOUtils.bytesToIntLittleEndian(bBuf, 4);
          FileUtils.skip(is, numChunks * 16L);
          pos += numChunks * 16L;
        }
        mLinearIndexPositions[seqNo] = pos;
        IOUtils.readFully(is, bBuf, 0, 4);
        pos += 4;
        final int numLinear = ByteArrayIOUtils.bytesToIntLittleEndian(bBuf, 0);
        FileUtils.skip(is, numLinear * 8L);
        pos += numLinear * 8L;
      }
    }
  }

  @Override
  public InputStream openIndexFile() throws IOException {
    return new FileInputStream(mIndexFile);
  }
}
