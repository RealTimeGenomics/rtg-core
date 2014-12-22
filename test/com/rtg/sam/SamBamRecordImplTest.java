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

import java.io.File;
import java.util.Arrays;

import com.rtg.util.io.FileUtils;
import com.rtg.util.test.FileHelper;

import junit.framework.TestCase;


/**
 * Uses the first SAM record set up by SamReaderTest to test the
 * extra named field getters of the SamBamRecord wrapper class.
 *
 */
public class SamBamRecordImplTest extends TestCase {


  public void testGetReadName() throws Exception {

    final File tmpDir = FileUtils.createTempDir("bamreader", "blah");
    try {
      final File samFile = new File(tmpDir, "sam.sam");
      final File bamFile = new File(tmpDir, "sam.bam");
      final File bamIndex = new File(tmpDir, "sam.bai");
      final String header = "@HD\tVN:1.0\tSO:coordinate\n"
        + "@SQ\tSN:simulatedSequence\tLN:10000\n"
        + "@PG\tID:rtg\tVN:v2.0-EAP2dev build 20721 (2009-10-01)\n";
      final String content = header
                       + "962\t163\tsimulatedSequence\t16\t255\t5H2I3=5P2M3X5D5N6S\t=\t133\t152\tGTTTCCTCNCCGTAGTGGAATCGATGCTAATGAGAC\taaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa\tAS:A:c\tNM:i:1\tMQ:f:2.55\tIH:i:1\tHH:Z:blah\n";
      FileUtils.stringToFile(content, samFile);
      Sam2Bam.convertSamToBam(bamFile, bamIndex, Arrays.asList(new File[] {samFile}));

      try (BamReader br = new BamReader(bamFile)) {
        assertTrue(br.hasNext());
        br.next();

        final SamBamRecordImpl sbri = new SamBamRecordImpl(br);

        assertNotNull(sbri.mReader);

        assertEquals("962", sbri.getReadName());
        assertEquals(962, sbri.getReadNumber());
        sbri.getFlags();
        assertTrue(sbri.hasAttribute("MQ"));
        assertEquals('f', sbri.getAttributeType("MQ"));
        assertEquals("blah", (String) sbri.getAttributeValue("HH"));
        assertEquals(5, sbri.getNumFields());
        assertEquals("5H2I3=5P2M3X5D5N6S", sbri.getField(SamBamConstants.CIGAR_FIELD));
        assertEquals(163, sbri.getIntField(1));

        assertEquals("[AS, NM, MQ, IH, HH]", Arrays.toString(sbri.getAttributeTags()));

        assertEquals(1, sbri.getIntAttribute("IH"));
        assertEquals(163, sbri.getFlags());

        try {
          sbri.getFieldNumFromTag("MQ");
          fail();
        } catch (final UnsupportedOperationException uoe) {
          assertEquals("Not supported.", uoe.getMessage());
        }
      }
    } finally {
      FileHelper.deleteAll(tmpDir);
    }
  }
}
