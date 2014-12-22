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

import java.io.IOException;
import java.util.Arrays;

import com.rtg.tabix.TabixIndexer.TabixOptions;
import com.rtg.util.io.ByteArrayIOUtils;
import com.rtg.util.io.IOUtils;

import net.sf.samtools.util.BlockCompressedInputStream;

/**
 * Helper class for header field encapsulation
 */
public class TabixHeader {

  private static final int FIXED_SIZE = 36;
  private final int mNumSequences;
  private final TabixIndexer.TabixOptions mOptions;
  private final byte[] mSequenceNames;
  private String[] mSequenceNamesUnpacked;

  TabixHeader(int numSequences, TabixIndexer.TabixOptions options, byte[] sequenceNames) {
    mNumSequences = numSequences;
    mOptions = options;
    mSequenceNames = sequenceNames;
    mSequenceNamesUnpacked = new String[mNumSequences];
    int seqNo = 0;
    int sp = 0;
    for (int i = 0; i < mSequenceNames.length; i++) {
      if (mSequenceNames[i] == 0) {
        mSequenceNamesUnpacked[seqNo++] = new String(mSequenceNames, sp, i - sp);
        sp = i + 1;
      }
    }
  }

  static TabixHeader readHeader(BlockCompressedInputStream is) throws IOException {
    final byte[] fixedData = new byte[FIXED_SIZE];
    IOUtils.readFully(is, fixedData, 0, FIXED_SIZE);
    final int numberReferences = ByteArrayIOUtils.bytesToIntLittleEndian(fixedData, 4);
    final int format = ByteArrayIOUtils.bytesToIntLittleEndian(fixedData, 8);
    final int seqCol = ByteArrayIOUtils.bytesToIntLittleEndian(fixedData, 12) - 1;
    final int begCol = ByteArrayIOUtils.bytesToIntLittleEndian(fixedData, 16) - 1;
    final int endCol = ByteArrayIOUtils.bytesToIntLittleEndian(fixedData, 20) - 1;
    final int meta = ByteArrayIOUtils.bytesToIntLittleEndian(fixedData, 24);
    final int skip = ByteArrayIOUtils.bytesToIntLittleEndian(fixedData, 28);
    final int sequenceNameLength = ByteArrayIOUtils.bytesToIntLittleEndian(fixedData, 32);
    final byte[] sequenceNames = new byte[sequenceNameLength];
    IOUtils.readFully(is, sequenceNames, 0, sequenceNameLength);
    return new TabixHeader(numberReferences, new TabixIndexer.TabixOptions(format, seqCol, begCol, endCol, meta, skip), sequenceNames);
  }

  static TabixHeader mergeHeaders(TabixHeader firstHeader, TabixHeader nextHeader) {
    int startPos = 0;
    int numSequences;
    if (firstHeader.mNumSequences > 0 && nextHeader.mNumSequences > 0
            && firstHeader.mSequenceNamesUnpacked[firstHeader.mNumSequences - 1].equals(nextHeader.mSequenceNamesUnpacked[0])) {
      final byte[] secondNames = nextHeader.mSequenceNames;
      for (int i = 0; i < secondNames.length; i++) {
        if (secondNames[i] == 0) {
          startPos = i + 1;
          break;
        }
      }
      numSequences = firstHeader.mNumSequences + nextHeader.mNumSequences - 1;
    } else {
      numSequences = firstHeader.mNumSequences + nextHeader.mNumSequences;
    }
    final byte[] mergedNames = Arrays.copyOf(firstHeader.mSequenceNames, firstHeader.mSequenceNames.length + nextHeader.mSequenceNames.length - startPos);
    System.arraycopy(nextHeader.mSequenceNames, startPos, mergedNames, firstHeader.mSequenceNames.length, nextHeader.mSequenceNames.length - startPos);
    return new TabixHeader(numSequences, firstHeader.mOptions, mergedNames);
  }

  public int getNumSequences() {
    return mNumSequences;
  }

  public TabixOptions getOptions() {
    return mOptions;
  }

  public String[] getSequenceNamesUnpacked() {
    return mSequenceNamesUnpacked;
  }

}
