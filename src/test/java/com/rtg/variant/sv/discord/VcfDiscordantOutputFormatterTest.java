/*
 * Copyright (c) 2018. Real Time Genomics Limited.
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the
 *    distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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
import com.rtg.variant.sv.bndeval.BreakpointGeometry;
import com.rtg.variant.sv.bndeval.Orientation;

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
