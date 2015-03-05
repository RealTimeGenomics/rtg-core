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

import java.io.File;
import java.io.IOException;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import htsjdk.samtools.util.BlockCompressedInputStream;

/**
 * Grab the record with the provided coordinates
 */
public final class TabixSniper {

  private TabixSniper() { }

  /**
   * @param args command line arguments
   * @throws IOException if an IO error occurs
   */
  public static void main(String[] args) throws IOException {
    final BlockCompressedLineReader bclr = new BlockCompressedLineReader(new BlockCompressedInputStream(new File(args[0])));
    final Pattern p = Pattern.compile("\\(([0-9]+), ([0-9]+)\\)");
    final Matcher m = p.matcher(args[1]);
    m.find();
    final long big = Long.parseLong(m.group(1));
    final long small = Long.parseLong(m.group(2));
    bclr.seek((big << 16L) + small);
    int num = args.length > 2 ? Integer.parseInt(args[2]) : 1;
    while (num > 0) {
      System.out.println(bclr.readLine());
      num--;
    }
  }
}
