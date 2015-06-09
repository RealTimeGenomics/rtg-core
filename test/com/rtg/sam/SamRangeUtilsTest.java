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

import static com.rtg.util.StringUtils.TAB;

import java.io.ByteArrayInputStream;
import java.io.File;
import java.io.IOException;

import com.rtg.util.diagnostic.NoTalkbackSlimException;
import com.rtg.util.intervals.RangeList;
import com.rtg.util.intervals.ReferenceRanges;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.TestDirectory;

import htsjdk.samtools.SAMFileHeader;
import junit.framework.TestCase;

/**
 * Tests for SamRangeUtils
 */
public class SamRangeUtilsTest extends TestCase {

  private static final String NL = "\n";  // SAM files always use \n.

  public SamRangeUtilsTest(final String name) {
    super(name);
  }

  private static final String SAM_HEADER = ""
    + "@HD" + TAB + "VN:1.0" + TAB + "SO:coordinate" + NL
    + "@SQ" + TAB + "SN:g0" + TAB + "LN:200" + NL
    + "@SQ" + TAB + "SN:g1" + TAB + "LN:400" + NL
    + "@SQ" + TAB + "SN:g2" + TAB + "LN:600" + NL
    ;

  private static final String BED_REGIONS = ""
    + "g1\t0\t50\n"
    + "g1\t100\t500\n"
    + "g2\t20\t500\n"
    ;

  public void testCreateSingleRangeList() throws IOException {
    final ByteArrayInputStream bis = new ByteArrayInputStream(SAM_HEADER.getBytes());
    final SAMFileHeader header = SamUtils.getSingleHeader(bis);

    ReferenceRanges<String> refRanges = SamRangeUtils.createExplicitReferenceRange(header, new SamRegionRestriction("g1"));

    assertNull(refRanges.get("g0"));
    assertNull(refRanges.get("g2"));
    assertNotNull(refRanges.get("g1"));

    RangeList<String> seqRanges = refRanges.get("g1");
    assertNotNull(seqRanges.find(100));
    assertNull(seqRanges.find(1000));

    try {
      SamRangeUtils.createExplicitReferenceRange(header, new SamRegionRestriction("g3"));
    } catch (NoTalkbackSlimException e) {
      //
    }
  }

  public void testCreateBedRangeList() throws IOException {
    try (final TestDirectory dir = new TestDirectory("samrangeutils")) {
      final ByteArrayInputStream bis = new ByteArrayInputStream(SAM_HEADER.getBytes());
      final SAMFileHeader header = SamUtils.getSingleHeader(bis);
      final File bedfile = new File(dir, "regions.bed");
      FileUtils.stringToFile(BED_REGIONS, bedfile);
      ReferenceRanges<String> refRanges = SamRangeUtils.createBedReferenceRanges(bedfile);
      assertNull(refRanges.get("g0"));
      assertNotNull(refRanges.get("g1"));
      assertNotNull(refRanges.get("g2"));

      RangeList<String> seqRanges = refRanges.get("g1");
      assertNotNull(seqRanges.find(100));
      assertNull(seqRanges.find(1000));
      assertNotNull(seqRanges.find(450));

      refRanges = SamRangeUtils.createBedReferenceRanges(header, bedfile);
      seqRanges = refRanges.get("g1");
      assertNull(seqRanges.find(450)); // Validation against header will have clipped to length of g1
    }
  }

  public void testCreateFullRangeList() throws IOException {
    final ByteArrayInputStream bis = new ByteArrayInputStream(SAM_HEADER.getBytes());
    final SAMFileHeader header = SamUtils.getSingleHeader(bis);

    ReferenceRanges<String> refRanges = SamRangeUtils.createFullReferenceRanges(header);

    assertNotNull(refRanges.get("g0"));
    assertNotNull(refRanges.get("g1"));
    assertNotNull(refRanges.get("g2"));
    assertNull(refRanges.get("g3"));

    RangeList<String> seqRanges = refRanges.get("g1");
    assertNotNull(seqRanges.find(100));
    assertNull(seqRanges.find(1000));
  }

  public void testConvertNameToIdKeys() throws IOException {
    final ByteArrayInputStream bis = new ByteArrayInputStream(SAM_HEADER.getBytes());
    final SAMFileHeader header = SamUtils.getSingleHeader(bis);

    ReferenceRanges<String> refRanges = SamRangeUtils.createFullReferenceRanges(header);
    refRanges.setIdMap(SamUtils.getSequenceIdLookup(header.getSequenceDictionary()));

    assertNotNull(refRanges.get(0));
    assertNotNull(refRanges.get(1));
    assertNotNull(refRanges.get(2));
    assertNull(refRanges.get(3));

    RangeList<String> seqRanges = refRanges.get(1);
    assertNotNull(seqRanges.find(100));
    assertNull(seqRanges.find(1000));
  }

  public void testGetReferenceRanges() throws IOException {
    final ByteArrayInputStream bis = new ByteArrayInputStream(SAM_HEADER.getBytes());
    final SAMFileHeader header = SamUtils.getSingleHeader(bis);

    SamFilterParams params = new SamFilterParams.SamFilterParamsBuilder().create();
    ReferenceRanges<String> refRanges = SamRangeUtils.createReferenceRanges(header, params);
    assertEquals(3, refRanges.sequenceNames().size());

  }

}
