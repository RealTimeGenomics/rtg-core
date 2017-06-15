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

package com.rtg.variant.sv.discord;

import static com.rtg.util.StringUtils.LS;
import static com.rtg.util.StringUtils.TAB;

import java.io.File;
import java.io.IOException;

import com.rtg.AbstractTest;
import com.rtg.reader.ReaderTestUtils;
import com.rtg.reader.SdfId;
import com.rtg.reader.SequencesReader;
import com.rtg.util.io.TestDirectory;

/**
 */
public class VcfDiscordantOutputFormatterTest extends AbstractTest {

  private static final String REF = ">f" + LS
  //1234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890
  + "actgactactactgaccccccccatgtactgacgtcttagctgatgacgtaccatactgacgtatgacgtactgactgacgtactactgactgactgaca"
  + "actgactactactgaccccccccatgtactgacgtcttagctgatgacgtaccatactgacgtatgacgtactgactgacgtactactgactgactgaca"
  + "actgactactactgaccccccccatgtactgacgtcttagctgatgacgtaccatactgacgtatgacgtactgactgacgtactactgactgactgaca" + LS
  + ">s" + LS
  //1234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890
  + "actgactactactgaccccccccatgtactgacgtcttagctgatgacgtaccatactgacgtatgacgtactgactgacgtactactgactgactgaca"
  + "actgactactactgaccccccccatgtactgacgtcttagctgatgacgtaccatactgacgtatgacgtactgactgacgtactactgactgactgaca"
  + "actgactactactgaccccccccatgtactgacgtcttagctgatgacgtaccatactgacgtatgacgtactgactgacgtactactgactgactgaca"
  + LS;

  public void test() throws IOException {
    try (final TestDirectory t = new TestDirectory("vcfdiscord")) {
      try (SequencesReader templateReader = ReaderTestUtils.getReaderDNA(REF, new File(t, "ref"), new SdfId())) {
        final BreakpointConstraint bg = new BreakpointConstraint(new BreakpointGeometry(Orientation.UU, "f", "s", 42, 50, 142, 150, 190, 198), 0, 10.0);
        final DiscordantReadSet drs = new DiscordantReadSet("f", 60, bg);
        assertEquals(42, drs.unionPosition());
        assertEquals(102, drs.flushPosition());
        final String exp = "f" + TAB
          + "45" + TAB
          + "." + TAB
          + "T" + TAB
          + "T]s:146]" + TAB
          + "." + TAB
          + "PASS" + TAB
          + "IMPRECISE;SVTYPE=BND;DP=1;CIPOS=-3,3;CV=0;AR=0.0" + TAB //TODO check CIPOS
          + "GT" + TAB
          + "1/1";
        assertEquals(exp, new VcfDiscordantOutputFormatter(templateReader).vcfRecord(drs, 0, 0.0).toString());
      }
    }
  }

  public void testNegative() throws IOException {
    try (final TestDirectory t = new TestDirectory("vcfdiscord")) {
      try (SequencesReader templateReader = ReaderTestUtils.getReaderDNA(REF, new File(t, "ref"), new SdfId())) {

        final BreakpointConstraint bg = new BreakpointConstraint(new BreakpointGeometry(Orientation.UD, "f", "f", -20, 2, 3, -40, -23, 42), 0, 10.0); //new BreakpointConstraint(new BreakpointGeometry(Orientation.UU, "f", "s", 42, 50, 142, 150, 190, 198), 0, 10.0);
        final DiscordantReadSet drs = new DiscordantReadSet("f", 60, bg);
        assertEquals(-20, drs.unionPosition());
        assertEquals(40, drs.flushPosition());
        final String exp = "f" + TAB
          + "1" + TAB
          + "." + TAB
          + "A" + TAB
          + "A[f:1[" + TAB
          + "." + TAB
          + "PASS" + TAB
          + "IMPRECISE;SVTYPE=BND;DP=1;CIPOS=-21,1;CV=0;AR=0.0" + TAB //TODO check CIPOS
          + "GT" + TAB
          + "1/1";
        assertEquals(exp, new VcfDiscordantOutputFormatter(templateReader).vcfRecord(drs, 0, 0.0).toString());
      }
    }
  }

  public void testInconsistent() throws IOException {
    try (final TestDirectory t = new TestDirectory("vcfdiscord")) {
      try (SequencesReader templateReader = ReaderTestUtils.getReaderDNA(REF, new File(t, "ref"), new SdfId())) {
        final BreakpointConstraint bg = new BreakpointConstraint(new BreakpointGeometry(Orientation.UU, "f", "s", 42, 50, 142, 150, 190, 198), 0, 10.0);
        final DiscordantReadSet drs = new DiscordantReadSet("f", 60, bg);
        final BreakpointConstraint bg2 = new BreakpointConstraint(new BreakpointGeometry(Orientation.UU, "f", "s", 842, 850, 942, 950, 1792, 1796), 550, 50.0);
        final DiscordantReadSet drs2 = new DiscordantReadSet("f", 860, bg2);
        drs.addAll(drs2);
        assertEquals(42, drs.unionPosition());
        assertEquals(102, drs.flushPosition());
        final String exp = "f" + TAB
          + "45" + TAB
          + "." + TAB
          + "T" + TAB
          + "T]s:146]" + TAB
          + "." + TAB
          + "INCONSISTENT" + TAB
          + "IMPRECISE;SVTYPE=BND;DP=2;CIPOS=-3,3;CV=0;AR=0.0" + TAB //TODO check CIPOS
          + "GT" + TAB
          + "1/1";
        assertEquals(exp, new VcfDiscordantOutputFormatter(templateReader).vcfRecord(drs, 0, 0.0).toString());
      }
    }
  }

}
