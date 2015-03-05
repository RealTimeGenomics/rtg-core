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

import java.io.BufferedInputStream;
import java.io.EOFException;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.Map;
import java.util.SortedMap;
import java.util.TreeMap;

import com.rtg.util.StringUtils;
import com.rtg.util.io.ByteArrayIOUtils;
import com.rtg.util.io.IOUtils;

import htsjdk.samtools.util.BlockCompressedInputStream;

/**
 * Useful for testing BAM/TABIX indices
 */
public final class IndexTestUtils {


  private IndexTestUtils() {

  }

  private enum TbiFormat {
    GENERIC, SAM, VCF
  }

  /**
   * Creates a string representation of the TABIX index
   * @param is stream
   * @return string
   * @throws IOException if an IO error occurs
   */
  public static String tbiIndexToUniqueString(InputStream is) throws IOException {
    final StringBuilder ret = new StringBuilder();
    final byte[] buf = new byte[4096];
    readIOFully(is, buf, 4);
    final String header = new String(buf, 0, 4);
    ret.append("Header correct: ").append(header.equals("TBI\u0001")).append(StringUtils.LS);
    readIOFully(is, buf, 4);
    final int numRefs = ByteArrayIOUtils.bytesToIntLittleEndian(buf, 0);
    ret.append("numRefs: ").append(numRefs).append(StringUtils.LS);
    readIOFully(is, buf, 28);
    final int format = ByteArrayIOUtils.bytesToIntLittleEndian(buf, 0);
    final int colSeq = ByteArrayIOUtils.bytesToIntLittleEndian(buf, 4);
    final int colBeg = ByteArrayIOUtils.bytesToIntLittleEndian(buf, 8);
    final int colEnd = ByteArrayIOUtils.bytesToIntLittleEndian(buf, 12);
    final int meta = ByteArrayIOUtils.bytesToIntLittleEndian(buf, 16);
    final int skip = ByteArrayIOUtils.bytesToIntLittleEndian(buf, 20);
    final int refNameLength = ByteArrayIOUtils.bytesToIntLittleEndian(buf, 24);
    final String formatStr;
    formatStr = TbiFormat.values()[format & 0xffff].name();
    ret.append("Format: ").append(formatStr).append(" 0-based: ").append((format & 0x10000) != 0).append(StringUtils.LS);
    ret.append("Columns: (refName:Start-End) ").append(colSeq).append(":").append(colBeg).append("-").append(colEnd).append(StringUtils.LS);
    ret.append("Meta: ").append((char) meta).append(StringUtils.LS);
    ret.append("Skip: ").append(skip).append(StringUtils.LS);
    final byte[] names = new byte[refNameLength];
    readIOFully(is, names, names.length);
    ret.append("Sequence names: ");
    boolean first = true;
    int off = 0;
    for (int i = 0; i < numRefs; i++) {
      int newOff = off;
      while (newOff < names.length && names[newOff] != 0) {
        newOff++;
      }
      if (!first) {
        ret.append(", ");
      }
      ret.append(new String(names, off, newOff - off));
      off = newOff + 1;
      first = false;
    }
    ret.append(StringUtils.LS);
    ret.append(indicesToUniqueString(is, numRefs)).append(StringUtils.LS);
    return ret.toString();
  }

  /**
   * Creates a string representation of the BAM index
   * @param is stream
   * @return string
   * @throws IOException if an IO error occurs
   */
  public static String bamIndexToUniqueString(InputStream is) throws IOException {
    final StringBuilder ret = new StringBuilder();
    final byte[] buf = new byte[4096];
    readIOFully(is, buf, 4);
    final String header = new String(buf, 0, 4);
    ret.append("Header correct: ").append(header.equals("BAI\u0001")).append(StringUtils.LS);
    readIOFully(is, buf, 4);
    final int numRefs = ByteArrayIOUtils.bytesToIntLittleEndian(buf, 0);
    ret.append("numRefs: ").append(numRefs).append(StringUtils.LS);
    ret.append(indicesToUniqueString(is, numRefs));
    final int len = is.read(buf, 0, 8);
    if (len > 0) {
      IOUtils.readFully(is, buf, len, 8 - len);
      ret.append("Number of unmapped records with no coordinates: ").append(ByteArrayIOUtils.bytesToLongLittleEndian(buf, 0));
      ret.append(StringUtils.LS);
    }
    return ret.toString();
  }

  private static String indicesToUniqueString(InputStream is, int numRefs) throws IOException {
    final StringBuilder ret = new StringBuilder();
    final byte[] buf = new byte[4096];
    for (int i = 0; i < numRefs; i++) {
      ret.append("refId: ").append(i).append(StringUtils.LS);
      readIOFully(is, buf, 4);
      final int numBins = ByteArrayIOUtils.bytesToIntLittleEndian(buf, 0);
      ret.append("  numBins: ").append(numBins).append(StringUtils.LS);
      //final SortedMap<Integer, ArrayList<LongPair>> binSet = new TreeMap<>();
      final SortedMap<Integer, VirtualOffsets> binSet = new TreeMap<>();
      for (int j = 0; j < numBins; j++) {
        readIOFully(is, buf, 8);
        final int bin = ByteArrayIOUtils.bytesToIntLittleEndian(buf, 0);
        final int numChunks = ByteArrayIOUtils.bytesToIntLittleEndian(buf, 4);
        final VirtualOffsets chunks = new VirtualOffsets();
        binSet.put(bin, chunks);
        for (int k = 0; k < numChunks; k++) {
          readIOFully(is, buf, 16);
          chunks.add(ByteArrayIOUtils.bytesToLongLittleEndian(buf, 0), ByteArrayIOUtils.bytesToLongLittleEndian(buf, 8), null);
        }
      }
      boolean start = false;
      int cc = 0;
      for (Map.Entry<Integer, VirtualOffsets> entry : binSet.entrySet()) {
        if (cc == 0) {
          if (start) {
            ret.append(StringUtils.LS);
          }
          start = true;
          ret.append("    binSet: ");
        }
        ret.append(entry).append(", ");
        cc++;
        if (cc == 10) {
          cc = 0;
        }
      }
      ret.append(StringUtils.LS);
      //+ binSet + StringUtils.LS);
      readIOFully(is, buf, 4);
      final int numIntervals = ByteArrayIOUtils.bytesToIntLittleEndian(buf, 0);
      final long[] intervals = new long[numIntervals];
      ret.append("  numIntervals: ").append(numIntervals).append(StringUtils.LS);
      for (int j = 0; j < numIntervals; j++) {
        readIOFully(is, buf, 8);
        final long off = ByteArrayIOUtils.bytesToLongLittleEndian(buf, 0);
        intervals[j] = off;
      }
      start = false;
      cc = 0;
      int index = 0;
      for (long off : intervals) {
        if (cc == 0) {
          if (start) {
            ret.append(StringUtils.LS);
          }
          start = true;
          ret.append("    intervals ").append(index).append("-").append(index + 9).append(": ");
        }
        ret.append(VirtualOffsets.offsetToString(off)).append(", ");
        cc++;
        if (cc == 10) {
          cc = 0;
        }
        index++;
      }
      ret.append(StringUtils.LS);
    }
    return ret.toString();
  }

  private static int readIOFully(InputStream is, byte[] buf, int amount) throws IOException {
    int tot = 0;
    if (amount != 0) {
      int len;
      do {
        len = is.read(buf, tot, amount - tot);
        if (len <= 0) {
          throw new EOFException();
        }
        tot += len;
      } while (tot < amount);
    }
    return tot;
  }

  /**
   * prints the index file in string form
   * @param args filename of index
   * @throws IOException if an IO error occurs.
   */
  public static void main(String[] args) throws IOException {
    InputStream is = new BufferedInputStream(new FileInputStream(args[0]));
    try {
      if (BlockCompressedInputStream.isValidFile(is)) {
        is = new BlockCompressedInputStream(is);
        System.out.println(tbiIndexToUniqueString(is));
      } else {
        System.out.println(bamIndexToUniqueString(is));
      }
    } finally {
      is.close();
    }
  }
}
