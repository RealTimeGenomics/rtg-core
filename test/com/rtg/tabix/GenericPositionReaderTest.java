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

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.IOException;

import com.rtg.util.StringUtils;

import net.sf.samtools.util.BlockCompressedInputStream;
import net.sf.samtools.util.BlockCompressedOutputStream;

import junit.framework.TestCase;

/**
 * Test class
 */
public class GenericPositionReaderTest extends TestCase {

 private static final String SNP_FILE = "#name\tposition\ttype\treference\tprediction\tposterior\tcoverage\tcorrection\tsupport_statistics" + StringUtils.LS
          + "simulatedSequence1\t184\te\tT\tA:T\t5.0\t19\t0.378\tA\t6\t0.119\tT\t13\t0.259" + StringUtils.LS
          + "simulatedSequence1\t2180\te\tC\tG:T\t28.7\t35\t0.697\tG\t18\t0.358\tT\t17\t0.338" + StringUtils.LS;
  public void testSomeMethod() throws IOException {
    final ByteArrayOutputStream baos = new ByteArrayOutputStream();
    final BlockCompressedOutputStream out = new BlockCompressedOutputStream(baos, null);
    out.write(SNP_FILE.getBytes());
    out.close();
    final ByteArrayInputStream bais = new ByteArrayInputStream(baos.toByteArray());
    try (GenericPositionReader gpr = new GenericPositionReader(new BlockCompressedLineReader(new BlockCompressedInputStream(bais)), new TabixIndexer.TabixOptions(TabixIndexer.TabixOptions.FORMAT_GENERIC, 0, 1, 1, '#', 0, false))) {
      final int[] pos = {183, 2179};
      final int[] len = {1, 1};
      int i = 0;
      while (gpr.hasNext()) {
        gpr.next();
        assertEquals("simulatedSequence1", gpr.getReferenceName());
        assertEquals(pos[i], gpr.getStartPosition());
        assertEquals(len[i], gpr.getLengthOnReference());
        i++;
      }
      assertEquals(2, i);
    }
  }

}
