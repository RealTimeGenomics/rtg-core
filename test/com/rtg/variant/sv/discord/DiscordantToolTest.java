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

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.Map;
import java.util.SortedSet;
import java.util.TreeSet;

import com.rtg.launcher.OutputParams;
import com.rtg.launcher.SequenceParams;
import com.rtg.reader.ReaderTestUtils;
import com.rtg.util.InvalidParamsException;
import com.rtg.util.StringUtils;
import com.rtg.util.TestUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.io.FileUtils;
import com.rtg.util.test.FileHelper;
import com.rtg.variant.sv.ReadGroupStats;
import com.rtg.variant.sv.discord.DiscordantToolParams.DiscordantToolParamsBuilder;

import junit.framework.TestCase;

/**
 */
public class DiscordantToolTest extends TestCase {
  private static final String TAB = "\t";

  protected DiscordantToolParams makeParams(File ref, File input, File output, int minDepth, boolean intersectionOnly, Double ambiguity, Integer coverage) throws InvalidParamsException {
    final DiscordantToolParamsBuilder builder = DiscordantToolParams.builder();
    builder.name("Foo");
    builder.debugOutput(true);
    builder.genome(SequenceParams.builder().directory(ref).loadNames(true).useMemReader(true).create());
    final boolean gzip = false;
    final OutputParams outParams = new OutputParams(output, false, gzip);
    builder.outputParams(outParams);
//    builder.maxAmbiguity(ambiguity);
//    builder.maxCoverage(coverage);

    final Map<String, ReadGroupStats> rgs = new HashMap<>();
    rgs.put("RG1", new ReadGroupStats("RG1", 6, 70, 5, 50, 5, 5, 1, 1, 1, 1));
    builder.readGroupStatistics(rgs);
    builder.minBreakpointDepth(minDepth);
    builder.intersectionOnly(intersectionOnly);

    builder.ioThreads(1);
    final Collection<File> inputFiles = new ArrayList<>();
    inputFiles.add(input);
    builder.mapped(inputFiles);
    return builder.create();
  }

  static final String SAM = ""
      + "@HD" + TAB + "VN:1.0" + TAB + "SO:coordinate" + LS
      + "@SQ" + TAB + "SN:g1" + TAB + "LN:20" + LS
      + "@SQ" + TAB + "SN:g2" + TAB + "LN:20" + LS
      + "@RG" + TAB + "ID:RG1" + TAB + "PL:ILLUMINA" + TAB + "SM:bar" + LS
      + "a0" + TAB + "97" + TAB + "g1" + TAB +  "4" + TAB + "255" + TAB + "3=2I13=" + TAB + "g2" + TAB + "14" + TAB + "1" + TAB + "cgtagagagagagatgct" + TAB + "``````````````````" + TAB + "AS:i:3" + TAB + "RG:Z:RG1" + LS
      + "a1" + TAB + "97" + TAB + "g1" + TAB +  "5" + TAB + "255" + TAB + "3=2I13=" + TAB + "g2" + TAB + "15" + TAB + "1" + TAB + "cgtagagagagagatgct" + TAB + "``````````````````" + TAB + "AS:i:3" + TAB + "RG:Z:RG1" + LS
      + "a2" + TAB + "97" + TAB + "g1" + TAB +  "6" + TAB + "255" + TAB + "3=2I13=" + TAB + "g1" + TAB + "16" + TAB + "1" + TAB + "cgtagagagagagatgct" + TAB + "``````````````````" + TAB + "AS:i:3" + TAB + "RG:Z:RG1" + LS
      + "b2" + TAB + "145" + TAB + "g1" + TAB +  "16" + TAB + "255" + TAB + "12=2I4=" + TAB + "g1" + TAB + "6" + TAB + "1" + TAB + "cgtagagagagagatgct" + TAB + "``````````````````" + TAB + "AS:i:3" + TAB + "RG:Z:RG1" + LS
      + "b0" + TAB + "145" + TAB + "g2" + TAB +  "14" + TAB + "255" + TAB + "12=2I4=" + TAB + "g1" + TAB + "4" + TAB + "1" + TAB + "cgtagagagagagatgct" + TAB + "``````````````````" + TAB + "AS:i:3" + TAB + "RG:Z:RG1" + LS
      + "b1" + TAB + "145" + TAB + "g2" + TAB +  "15" + TAB + "255" + TAB + "12=2I4=" + TAB + "g1" + TAB + "5" + TAB + "1" + TAB + "cgtagagagagagatgct" + TAB + "``````````````````" + TAB + "AS:i:3" + TAB + "RG:Z:RG1" + LS
      ;

  static final String REF = ""
      //12345678901234567890
      + ">g1" + LS + "catgcagtcatactgactgtcatgcagtcatactgactgtcatgcagtcatactgactgtgctgtgtttgttgggtgggtttcgtcaaacacccacaaaatgtgttgtgacgtagatattagatagaagagtat" + LS
      + ">g2" + LS + "catgcagtcatactgactgtcatgcagtcatactgactgtcatgcagtcatactgactgtgctgtgtttgttgggtgggtttcgtcaaacacccacaaaatgtgttgtgacgtagatattagatagaagagtat" + LS
      + ">g3" + LS + "catgcagtcatactgactgtcatgcagtcatactgactgtcatgcagtcatactgactgtgctgtgtttgttgggtgggtttcgtcaaacacccacaaaatgtgttgtgacgtagatattagatagaagagtat" + LS;

  static final String EXP_OUT = ""
      + "2 g1 17 88 g2 14 -57 34 74 g1 18 87 g2 13 -56 34 74" + LS
      + "1 g1 18 -52 g1 22 92 34 74 g1 18 -52 g1 22 92 34 74" + LS
      + "1 g1 19 89 g1 15 -55 34 74 g1 19 89 g1 15 -55 34 74" + LS
      + "2 g2 17 -54 g1 20 91 34 74 g2 16 -53 g1 21 90 34 74"
      ;
  public void endToEnd(int minDepth, String expected, String sam, boolean intersectionOnly) throws Exception {
    Diagnostic.setLogStream();
    final File dir = FileHelper.createTempDirectory();
    try {
      final File seq = new File(dir, "ref.sdf");
      ReaderTestUtils.getDNADir(REF, seq);
      final File input = new File(dir, "input");
      FileUtils.stringToFile(sam, input);
      final File output = new File(dir, "output");
      final DiscordantToolParams p = makeParams(seq, input, output, minDepth, intersectionOnly, null, null);
      final DiscordantTool dt = new DiscordantTool(p, TestUtils.getNullOutputStream());
      dt.exec();
      assertTrue(output.exists());
      final File simple = new File(output, DiscordantTool.DEBUG_FILENAME);
      assertTrue(simple.exists());

      final String noheader = StringUtils.grep(FileUtils.fileToString(simple), "^[^#]").trim().replace(TAB, " ");
      //final String vcf = StringUtils.grep(FileUtils.fileToString(new File(output, DiscordantTool.OUTPUT_FILENAME)), "^[^#]").trim().replace(TAB, " ");
      //System.err.println(vcf);
      assertEquals(expected, noheader);
    } finally {
      FileHelper.deleteAll(dir);
    }
  }
  public void testEndToEnd() throws Exception {
    endToEnd(1, EXP_OUT, SAM, false);
  }

  public void testDepthCutoff() throws Exception {
    endToEnd(2, StringUtils.grep(EXP_OUT, "^2").trim(), SAM, false);
  }

  static final String SAM_HEADER = ""
      + "@HD" + TAB + "VN:1.0" + TAB + "SO:coordinate" + LS
      + "@SQ" + TAB + "SN:g1" + TAB + "LN:133" + LS
      + "@SQ" + TAB + "SN:g2" + TAB + "LN:133" + LS
      + "@SQ" + TAB + "SN:g3" + TAB + "LN:133" + LS
      + "@RG" + TAB + "ID:RG1" + TAB + "PL:ILLUMINA" + TAB + "SM:bar" + LS;
  static final String NULL_INTERSECT_SAM = ""
      + SAM_HEADER
      + "a0" + TAB + "97" + TAB + "g1" + TAB +  "4" + TAB + "255" + TAB + "9=" + TAB + "g2" + TAB + "14" + TAB + "1" + TAB + "cgtagagag" + TAB + "`````````" + TAB + "AS:i:0" + TAB + "RG:Z:RG1" + LS
      + "a1" + TAB + "97" + TAB + "g1" + TAB +  "5" + TAB + "255" + TAB + "3=2I13=" + TAB + "g2" + TAB + "15" + TAB + "1" + TAB + "cgtagagagagagatgct" + TAB + "``````````````````" + TAB + "AS:i:3" + TAB + "RG:Z:RG1" + LS
      + "a2" + TAB + "97" + TAB + "g1" + TAB +  "6" + TAB + "255" + TAB + "3=2I13=" + TAB + "g1" + TAB + "16" + TAB + "1" + TAB + "cgtagagagagagatgct" + TAB + "``````````````````" + TAB + "AS:i:3" + TAB + "RG:Z:RG1" + LS
//      + "a0b" + TAB + "97" + TAB + "g1" + TAB +  "30" + TAB + "255" + TAB + "9=" + TAB + "g2" + TAB + "44" + TAB + "1" + TAB + "cgtagagag" + TAB + "`````````" + TAB + "AS:i:0" + TAB + "RG:Z:RG1" + LS
      + "a3" + TAB + "97" + TAB + "g1" + TAB +  "50" + TAB + "255" + TAB + "9=" + TAB + "g2" + TAB + "64" + TAB + "1" + TAB + "cgtagagag" + TAB + "`````````" + TAB + "AS:i:0" + TAB + "RG:Z:RG1" + LS
      + "a4" + TAB + "97" + TAB + "g1" + TAB +  "60" + TAB + "255" + TAB + "9=" + TAB + "g2" + TAB + "74" + TAB + "1" + TAB + "cgtagagag" + TAB + "`````````" + TAB + "AS:i:0" + TAB + "RG:Z:RG1" + LS
      + "a5" + TAB + "97" + TAB + "g1" + TAB +  "70" + TAB + "255" + TAB + "9=" + TAB + "g2" + TAB + "84" + TAB + "1" + TAB + "cgtagagag" + TAB + "`````````" + TAB + "AS:i:0" + TAB + "RG:Z:RG1" + LS
      + "a6" + TAB + "97" + TAB + "g1" + TAB +  "70" + TAB + "255" + TAB + "9=" + TAB + "g3" + TAB + "84" + TAB + "1" + TAB + "cgtagagag" + TAB + "`````````" + TAB + "AS:i:0" + TAB + "RG:Z:RG1" + LS
      ;

  static final String NULL_INTERSECT_OUT = ""
      + "5 g1 13 149 g2 83 -57 26 74" + LS
      + "1 g1 19 89 g1 15 -55 34 74 g1 19 89 g1 15 -55 34 74" + LS // a2
      + "1 g1 79 149 g3 83 13 26 66 g1 79 149 g3 83 13 26 66" //a6
      ;
  public void testNullIntersection() throws Exception {
    endToEnd(0, NULL_INTERSECT_OUT, NULL_INTERSECT_SAM, false);
    endToEnd(0,  StringUtils.grepMinusV(NULL_INTERSECT_OUT, "g2").trim(), NULL_INTERSECT_SAM.trim(), true);
  }

  static final String FLUSH_SAM = ""
      + "@HD" + TAB + "VN:1.0" + TAB + "SO:coordinate" + LS
      + "@SQ" + TAB + "SN:g1" + TAB + "LN:100" + LS
      + "@SQ" + TAB + "SN:g2" + TAB + "LN:100" + LS
      + "@RG" + TAB + "ID:RG1" + TAB + "PL:ILLUMINA" + TAB + "SM:bar" + LS
      + "a0" + TAB + "97" + TAB + "g1" + TAB +  "4" + TAB + "255" + TAB + "3=2I13=" + TAB + "g2" + TAB + "14" + TAB + "1" + TAB + "cgtagagagagagatgct" + TAB + "``````````````````" + TAB + "AS:i:3" + TAB + "RG:Z:RG1" + LS
      + "a1" + TAB + "97" + TAB + "g1" + TAB +  "5" + TAB + "255" + TAB + "3=2I13=" + TAB + "g2" + TAB + "15" + TAB + "1" + TAB + "cgtagagagagagatgct" + TAB + "``````````````````" + TAB + "AS:i:3" + TAB + "RG:Z:RG1" + LS
      + "a2" + TAB + "97" + TAB + "g1" + TAB +  "6" + TAB + "255" + TAB + "3=2I13=" + TAB + "g1" + TAB + "16" + TAB + "1" + TAB + "cgtagagagagagatgct" + TAB + "``````````````````" + TAB + "AS:i:3" + TAB + "RG:Z:RG1" + LS
      + "b2" + TAB + "145" + TAB + "g1" + TAB +  "16" + TAB + "255" + TAB + "12=2I4=" + TAB + "g1" + TAB + "6" + TAB + "1" + TAB + "cgtagagagagagatgct" + TAB + "``````````````````" + TAB + "AS:i:3" + TAB + "RG:Z:RG1" + LS
      + "a4" + TAB + "97" + TAB + "g1" + TAB +  "70" + TAB + "255" + TAB + "3=2I13=" + TAB + "g1" + TAB + "16" + TAB + "1" + TAB + "cgtagagagagagatgct" + TAB + "``````````````````" + TAB + "AS:i:3" + TAB + "RG:Z:RG1" + LS
      + "b0" + TAB + "145" + TAB + "g2" + TAB +  "14" + TAB + "255" + TAB + "12=2I4=" + TAB + "g1" + TAB + "4" + TAB + "1" + TAB + "cgtagagagagagatgct" + TAB + "``````````````````" + TAB + "AS:i:3" + TAB + "RG:Z:RG1" + LS
      + "b1" + TAB + "145" + TAB + "g2" + TAB +  "15" + TAB + "255" + TAB + "12=2I4=" + TAB + "g1" + TAB + "5" + TAB + "1" + TAB + "cgtagagagagagatgct" + TAB + "``````````````````" + TAB + "AS:i:3" + TAB + "RG:Z:RG1" + LS
      ;

  static class FlushCheck extends DiscordantTool {
    StringBuilder mFlushPositions = new StringBuilder();
    FlushCheck(DiscordantToolParams p, OutputStream os) throws IOException {
      super(p, TestUtils.getNullOutputStream());
    }
    @Override
    public int flush(int start, int end) throws IOException {
      final int before = mReadSets.size();
      final int res = super.flush(start, end);
      mFlushPositions.append("Flushing: ").append(mTemplateName).append("[").append(start).append(":").append(end).append("], sizeBefore=").append(before).append(" sizeAfter=").append(mReadSets.size());
      mFlushPositions.append(StringUtils.LS);
      return res;
    }
  }

  static final String REF_FLUSH = ""
      //12345678901234567890
      + ">g1" + LS + "catgcagtcatactgactgtcatgcagtcatactgactgtcatgcagtcatactgactgtcatgcagtcatactgactgtcatgcagtcatactgactgtcatactgactgtcatgcagtcatactgactgt" + LS
      + ">g2" + LS + "catgcagtcatactgactgtcatgcagtcatactgactgtcatgcagtcatactgactgtcatgcagtcatactgactgtcatgcagtcatactgactgtcatactgactgtcatgcagtcatactgactgt" + LS;

  public void testFlushing() throws Exception {
    final File dir = FileHelper.createTempDirectory();
    try {
      final File seq = new File(dir, "ref.sdf");
      ReaderTestUtils.getDNADir(REF_FLUSH, seq);
      final File input = new File(dir, "input");
      FileUtils.stringToFile(FLUSH_SAM, input);
      final File output = new File(dir, "output");
      final DiscordantToolParams p = makeParams(seq, input, output, 1, false, null, null);
      final FlushCheck dt = new FlushCheck(p, TestUtils.getNullOutputStream());
      dt.exec();
      assertTrue(output.exists());
      final File simple = new File(output, DiscordantTool.OUTPUT_FILENAME);
      assertTrue(simple.exists());
      //System.err.println(dt.mFlushPositions.toString());
      TestUtils.containsAll(dt.mFlushPositions.toString()
          , "Flushing: g1[-1:-1], sizeBefore=0 sizeAfter=0"
          , "Flushing: g1[0:70], sizeBefore=4 sizeAfter=4"
          , "Flushing: g2[-1:132], sizeBefore=1 sizeAfter=0"
          , "Flushing: g2[0:-1], sizeBefore=0 sizeAfter=0"
          );
    } finally {
      FileHelper.deleteAll(dir);
    }
  }

  //single constraint
  public void testProcessConstraint0() {
    final BreakpointConstraint constraint = new BreakpointConstraint(new BreakpointGeometry(Orientation.UU, "x", "y", 10, 20, 30, 40, 42, 52), 41.0, 12.3);
    //System.err.println(constraint.gnuPlot());
    final SortedSet<DiscordantReadSet> readSets = new TreeSet<>(new DiscordantReadSet.FlushPositionComparator());
    DiscordantTool.processConstraint(constraint, readSets, "y", 100, null);
    final String exp = ""
        + "[DiscordantReadSet:" + LS
        + "union=Break-point constraint:UU x=10,20:x y=30,40:y r=42,52 Gap: mean=41.0 std.dev.=12.3" + LS
        + "intersection=Break-point constraint:UU x=10,20:x y=30,40:y r=42,52 Gap: mean=41.0 std.dev.=12.3" + LS
        + "constraints count=1" + LS
        + "    Break-point constraint:UU x=10,20:x y=30,40:y r=42,52 Gap: mean=41.0 std.dev.=12.3" + LS
        + "]"
        ;
    assertEquals(exp, readSets.toString());
  }

  //two disjoint constraints
  public void testProcessConstraint1() {
    final BreakpointConstraint c0 = new BreakpointConstraint(new BreakpointGeometry(Orientation.UU, "x", "y", 10, 20, 30, 40, 42, 52), 41.0, 12.3);
    final BreakpointConstraint c1 = new BreakpointConstraint(new BreakpointGeometry(Orientation.UU, "x", "y", 20, 30, 40, 50, 62, 72), 41.0, 12.3);
    //System.err.println(c0.gnuPlot());
    final SortedSet<DiscordantReadSet> readSets = new TreeSet<>(new DiscordantReadSet.FlushPositionComparator());
    DiscordantTool.processConstraint(c0, readSets, "x", 100, null);
    DiscordantTool.processConstraint(c1, readSets, "x", 100, null);
    final String exp = ""
        + "[DiscordantReadSet:" + LS
        + "union=Break-point constraint:UU x=10,20:x y=30,40:y r=42,52 Gap: mean=41.0 std.dev.=12.3" + LS
        + "intersection=Break-point constraint:UU x=10,20:x y=30,40:y r=42,52 Gap: mean=41.0 std.dev.=12.3" + LS
        + "constraints count=1" + LS
        + "    Break-point constraint:UU x=10,20:x y=30,40:y r=42,52 Gap: mean=41.0 std.dev.=12.3" + LS
        + ", DiscordantReadSet:" + LS
        + "union=Break-point constraint:UU x=20,30:x y=40,50:y r=62,72 Gap: mean=41.0 std.dev.=12.3" + LS
        + "intersection=Break-point constraint:UU x=20,30:x y=40,50:y r=62,72 Gap: mean=41.0 std.dev.=12.3" + LS
        + "constraints count=1" + LS
        + "    Break-point constraint:UU x=20,30:x y=40,50:y r=62,72 Gap: mean=41.0 std.dev.=12.3" + LS
        + "]"
        ;
    assertEquals(exp, readSets.toString());
  }

  //two overlapping constraints
  public void testProcessConstraint2() {
    final BreakpointConstraint c0 = new BreakpointConstraint(new BreakpointGeometry(Orientation.UU, "x", "y", 10, 20, 30, 40, 42, 52), 41.0, 12.3);
    final BreakpointConstraint c2 = new BreakpointConstraint(new BreakpointGeometry(Orientation.UU, "x", "y", 14, 24, 34, 44, 50, 60), 41.0, 12.3);
    //System.err.println(c0.gnuPlot());
    //System.err.println(c2.gnuPlot());
    final SortedSet<DiscordantReadSet> readSets = new TreeSet<>(new DiscordantReadSet.FlushPositionComparator());
    DiscordantTool.processConstraint(c0, readSets, "x", 100, null);
    DiscordantTool.processConstraint(c2, readSets, "x", 100, null);
    final String exp = ""
        + "[DiscordantReadSet:" + LS
        + "union=Break-point constraint:UU x=10,24:x y=30,44:y r=42,60 Gap: mean=41.0 std.dev.=8.7" + LS
        + "intersection=Break-point constraint:UU x=14,18:x y=34,38:y r=50,52 Gap: mean=41.0 std.dev.=8.7" + LS
        + "constraints count=2" + LS
        + "    Break-point constraint:UU x=10,20:x y=30,40:y r=42,52 Gap: mean=41.0 std.dev.=12.3" + LS
        + "    Break-point constraint:UU x=14,24:x y=34,44:y r=50,60 Gap: mean=41.0 std.dev.=12.3" + LS
        + "]"
        ;
    assertEquals(exp, readSets.toString());
  }

  //three overlapping constraints
  public void testProcessConstraint3() {
    final BreakpointConstraint c0 = new BreakpointConstraint(new BreakpointGeometry(Orientation.UU, "x", "y", 10, 20, 30, 40, 42, 52), 41.0, 12.3);
    final BreakpointConstraint c1 = new BreakpointConstraint(new BreakpointGeometry(Orientation.UU, "x", "y", 18, 28, 38, 48, 58, 68), 41.0, 12.3);
    final BreakpointConstraint c2 = new BreakpointConstraint(new BreakpointGeometry(Orientation.UU, "x", "y", 14, 24, 34, 44, 50, 60), 41.0, 12.3);
    //System.err.println(c0.gnuPlot());
    //System.err.println(c1.gnuPlot());
    //System.err.println(c2.gnuPlot());
    final SortedSet<DiscordantReadSet> readSets = new TreeSet<>(new DiscordantReadSet.FlushPositionComparator());
    DiscordantTool.processConstraint(c0, readSets, "x", 100, null);
    DiscordantTool.processConstraint(c1, readSets, "x", 100, null);
    DiscordantTool.processConstraint(c2, readSets, "x", 100, null);
    final String exp = ""
        + "[DiscordantReadSet:" + LS
        + "union=Break-point constraint:UU x=10,28:x y=30,48:y r=42,68 Gap: mean=41.0 std.dev.=7.1" + LS
        + "intersection=null" + LS
        + "constraints count=3" + LS
        + "    Break-point constraint:UU x=14,24:x y=34,44:y r=50,60 Gap: mean=41.0 std.dev.=12.3" + LS
        + "    Break-point constraint:UU x=10,20:x y=30,40:y r=42,52 Gap: mean=41.0 std.dev.=12.3" + LS
        + "    Break-point constraint:UU x=18,28:x y=38,48:y r=58,68 Gap: mean=41.0 std.dev.=12.3" + LS
        + "]"      ;
    assertEquals(exp, readSets.toString());
  }

}
