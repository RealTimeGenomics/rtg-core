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

import com.rtg.launcher.AbstractNanoTest;
import com.rtg.launcher.OutputParams;
import com.rtg.launcher.SequenceParams;
import com.rtg.reader.ReaderTestUtils;
import com.rtg.util.InvalidParamsException;
import com.rtg.util.StringUtils;
import com.rtg.util.TestUtils;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.TestDirectory;
import com.rtg.variant.sv.ReadGroupStats;
import com.rtg.variant.sv.discord.DiscordantToolParams.DiscordantToolParamsBuilder;

/**
 */
public class DiscordantToolTest extends AbstractNanoTest {
  private static final String TAB = "\t";

  protected DiscordantToolParams makeParams(File ref, File input, File output, int minDepth, boolean intersectionOnly) throws InvalidParamsException {
    final DiscordantToolParamsBuilder builder = DiscordantToolParams.builder();
    builder.name("Foo");
    builder.debugOutput(true);
    builder.genome(SequenceParams.builder().directory(ref).loadNames(true).useMemReader(true).create().readerParams());
    final boolean gzip = false;
    final OutputParams outParams = new OutputParams(output, false, gzip);
    builder.outputParams(outParams);
    builder.overlapFraction(0);
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
      + "@HD VN:1.0 SO:coordinate" + LS
      + "@SQ SN:g1 LN:20" + LS
      + "@SQ SN:g2 LN:20" + LS
      + "@RG ID:RG1 PL:ILLUMINA SM:bar" + LS
      + "a0 97 g1 4 255 3=2I13= g2 14 1 cgtagagagagagatgct `````````````````` AS:i:3 RG:Z:RG1" + LS
      + "a1 97 g1 5 255 3=2I13= g2 15 1 cgtagagagagagatgct `````````````````` AS:i:3 RG:Z:RG1" + LS
      + "a2 97 g1 6 255 3=2I13= g1 16 1 cgtagagagagagatgct `````````````````` AS:i:3 RG:Z:RG1" + LS
      + "b2 145 g1 16 255 12=2I4= g1 6 1 cgtagagagagagatgct `````````````````` AS:i:3 RG:Z:RG1" + LS
      + "b0 145 g2 14 255 12=2I4= g1 4 1 cgtagagagagagatgct `````````````````` AS:i:3 RG:Z:RG1" + LS
      + "b1 145 g2 15 255 12=2I4= g1 5 1 cgtagagagagagatgct `````````````````` AS:i:3 RG:Z:RG1" + LS
      ;

  static final String REF = ""
      //12345678901234567890
      + ">g1" + LS + "catgcagtcatactgactgtcatgcagtcatactgactgtcatgcagtcatactgactgtgctgtgtttgttgggtgggtttcgtcaaacacccacaaaatgtgttgtgacgtagatattagatagaagagtat" + LS
      + ">g2" + LS + "catgcagtcatactgactgtcatgcagtcatactgactgtcatgcagtcatactgactgtgctgtgtttgttgggtgggtttcgtcaaacacccacaaaatgtgttgtgacgtagatattagatagaagagtat" + LS
      + ">g3" + LS + "catgcagtcatactgactgtcatgcagtcatactgactgtcatgcagtcatactgactgtgctgtgtttgttgggtgggtttcgtcaaacacccacaaaatgtgttgtgacgtagatattagatagaagagtat" + LS;

  public void endToEnd(String label, int minDepth, String sam, boolean intersectionOnly) throws Exception {
    try (final TestDirectory dir = new TestDirectory()) {
      final File seq = new File(dir, "ref.sdf");
      ReaderTestUtils.getDNADir(REF, seq);
      final File input = new File(dir, "input.sam");
      FileUtils.stringToFile(sam.replaceAll(" ", TAB), input);
      final File output = new File(dir, "output");
      final DiscordantToolParams p = makeParams(seq, input, output, minDepth, intersectionOnly);
      final DiscordantTool dt = new DiscordantTool(p, TestUtils.getNullOutputStream());
      dt.exec();
      assertTrue(output.exists());
      final File simple = new File(output, DiscordantTool.DEBUG_FILENAME);
      assertTrue(simple.exists());

      final String noheader = StringUtils.grep(FileUtils.fileToString(simple), "^[^#]").trim().replace(TAB, " ");
      mNano.check(label, noheader);
    }
  }
  public void testEndToEnd() throws Exception {
    endToEnd("end-to-end.txt", 1, SAM, false);
  }

  public void testDepthCutoff() throws Exception {
    endToEnd("depth-2.txt", 2, SAM, false);
  }

  static final String SAM_HEADER = ""
      + "@HD VN:1.0 SO:coordinate" + LS
      + "@SQ SN:g1 LN:133" + LS
      + "@SQ SN:g2 LN:133" + LS
      + "@SQ SN:g3 LN:133" + LS
      + "@RG ID:RG1 PL:ILLUMINA SM:bar" + LS;
  static final String NULL_INTERSECT_SAM = ""
      + SAM_HEADER
      + "a0 97 g1 4 255 9= g2 14 1 cgtagagag ````````` AS:i:0 RG:Z:RG1" + LS
      + "a1 97 g1 5 255 3=2I13= g2 15 1 cgtagagagagagatgct `````````````````` AS:i:3 RG:Z:RG1" + LS
      + "a2 97 g1 6 255 3=2I13= g1 16 1 cgtagagagagagatgct `````````````````` AS:i:3 RG:Z:RG1" + LS
//      + "a0b 97 g1 30 255 9= g2 44 1 cgtagagag ````````` AS:i:0 RG:Z:RG1" + LS
      + "a3 97 g1 50 255 9= g2 64 1 cgtagagag ````````` AS:i:0 RG:Z:RG1" + LS
      + "a4 97 g1 60 255 9= g2 74 1 cgtagagag ````````` AS:i:0 RG:Z:RG1" + LS
      + "a5 97 g1 70 255 9= g2 84 1 cgtagagag ````````` AS:i:0 RG:Z:RG1" + LS
      + "a6 97 g1 70 255 9= g3 84 1 cgtagagag ````````` AS:i:0 RG:Z:RG1" + LS
      ;

  public void testNullIntersection() throws Exception {
    endToEnd("null-intersection-false.txt", 0, NULL_INTERSECT_SAM, false);
    endToEnd("null-intersection-true.txt", 0, NULL_INTERSECT_SAM.trim(), true);
  }

  static final String FLUSH_SAM = ""
      + "@HD VN:1.0 SO:coordinate" + LS
      + "@SQ SN:g1 LN:100" + LS
      + "@SQ SN:g2 LN:100" + LS
      + "@RG ID:RG1 PL:ILLUMINA SM:bar" + LS
      + "a0 97 g1 4 255 3=2I13= g2 14 1 cgtagagagagagatgct `````````````````` AS:i:3 RG:Z:RG1" + LS
      + "a1 97 g1 5 255 3=2I13= g2 15 1 cgtagagagagagatgct `````````````````` AS:i:3 RG:Z:RG1" + LS
      + "a2 97 g1 6 255 3=2I13= g1 16 1 cgtagagagagagatgct `````````````````` AS:i:3 RG:Z:RG1" + LS
      + "b2 145 g1 16 255 12=2I4= g1 6 1 cgtagagagagagatgct `````````````````` AS:i:3 RG:Z:RG1" + LS
      + "a4 97 g1 70 255 3=2I13= g1 16 1 cgtagagagagagatgct `````````````````` AS:i:3 RG:Z:RG1" + LS
      + "b0 145 g2 14 255 12=2I4= g1 4 1 cgtagagagagagatgct `````````````````` AS:i:3 RG:Z:RG1" + LS
      + "b1 145 g2 15 255 12=2I4= g1 5 1 cgtagagagagagatgct `````````````````` AS:i:3 RG:Z:RG1" + LS
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
      mFlushPositions.append(LS);
      return res;
    }
  }

  static final String REF_FLUSH = ""
      //12345678901234567890
      + ">g1" + LS + "catgcagtcatactgactgtcatgcagtcatactgactgtcatgcagtcatactgactgtcatgcagtcatactgactgtcatgcagtcatactgactgtcatactgactgtcatgcagtcatactgactgt" + LS
      + ">g2" + LS + "catgcagtcatactgactgtcatgcagtcatactgactgtcatgcagtcatactgactgtcatgcagtcatactgactgtcatgcagtcatactgactgtcatactgactgtcatgcagtcatactgactgt" + LS;

  public void testFlushing() throws Exception {
    try (final TestDirectory dir = new TestDirectory()) {
      final File seq = new File(dir, "ref.sdf");
      ReaderTestUtils.getDNADir(REF_FLUSH, seq);
      final File input = new File(dir, "input.sam");
      FileUtils.stringToFile(FLUSH_SAM.replaceAll(" ", TAB), input);
      final File output = new File(dir, "output");
      final DiscordantToolParams p = makeParams(seq, input, output, 1, false);
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
