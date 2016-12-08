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

package com.rtg.assembler;

import static com.rtg.assembler.PacBio.ERROR_FLOOR;
import static com.rtg.util.StringUtils.LS;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.rtg.assembler.graph.Graph;
import com.rtg.assembler.graph.MutableGraph;
import com.rtg.assembler.graph.implementation.GraphKmerAttribute;
import com.rtg.assembler.graph.implementation.PathArray;
import com.rtg.mode.DnaUtils;
import com.rtg.reader.ReaderTestUtils;
import com.rtg.util.TestUtils;
import com.rtg.util.io.MemoryPrintStream;

import junit.framework.TestCase;

/**
 */
public class PacBioTest extends TestCase {
  public void testEmptyHits() {
    final List<List<ContigPosition>> lists = new ArrayList<>();
    final PacBio.HitMap hitMap = PacBio.joinHits(lists, 2);
    assertEquals(0, hitMap.size());
  }
  public void testHitJoining() {
    final long[][] contigs = {{}, {}, {4, 1, 4, 5}, {}, {}, {}, {4, 5}, {1, 4}};
    final int[][] positions = {{}, {}, {4, 10, 19, 6}, {}, {}, {}, {9 - ERROR_FLOOR, 9 + ERROR_FLOOR}, {15 - ERROR_FLOOR, 24 + ERROR_FLOOR}};
    final PacBio.HitMap hitMap = mapHits(contigs, positions);
    assertEquals(3, hitMap.size());
    assertEquals("[HitCollection[Hit{C10, R2}], HitCollection[Hit{C" + (15 - ERROR_FLOOR) + ", R7}]]", hitMap.get(1L).toString());
    assertEquals("[HitCollection[Hit{C4, R2}, Hit{C" + (9 - ERROR_FLOOR) + ", R6}], HitCollection[Hit{C19, R2}], HitCollection[Hit{C" + (24 + ERROR_FLOOR) + ", R7}]]", hitMap.get(4L).toString());
    assertEquals("[HitCollection[Hit{C6, R2}, Hit{C" + (9 + ERROR_FLOOR) + ", R6}]]", hitMap.get(5L).toString());
  }

  public void testEqualOffset() {
    final long[][] contigs = {{}, {1}, {}, {1}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {1}};
    final int[][] positions = {{}, {3}, {}, {1}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {15}};
    final PacBio.HitMap hitMap = mapHits(contigs, positions);
    assertEquals(1, hitMap.size());
    assertEquals("[HitCollection[Hit{C3, R1}], HitCollection[Hit{C1, R3}], HitCollection[Hit{C15, R15}]]", hitMap.get(1L).toString());


  }
  public void testInsertOffset() {
    final long[][] contigs = {{}, {1}, {}, {1}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {1}, {}};
    final int[][] positions = {{}, {3}, {}, {1}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {15}, {}};
    final PacBio.HitMap hitMap = mapHits(contigs, positions);
    assertEquals(1, hitMap.size());
    assertEquals("[HitCollection[Hit{C3, R1}, Hit{C15, R14}], HitCollection[Hit{C1, R3}]]", hitMap.get(1L).toString());
  }
  public void testDeleteOffset() {
    final long[][] contigs = {{}, {1}, {}, {1}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {1}};
    final int[][] positions = {{}, {3}, {}, {1}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {14}};
    final PacBio.HitMap hitMap = mapHits(contigs, positions);
    assertEquals(1, hitMap.size());
    assertEquals("[HitCollection[Hit{C3, R1}], HitCollection[Hit{C1, R3}, Hit{C14, R15}]]", hitMap.get(1L).toString());
  }

  public void testRepeatBlock() {
    final long[][] contigs = {{}, {1, 1, 1, 1}, {1, 1, 1, 1}, {1, 1, 1, 1}, {1, 1, 1, 1}, {1, 1, 1, 1}, {1, 1, 1, 1}, {1, 1, 1, 1}, {1, 1, 1, 1}, {1, 1, 1, 1}, {1, 1, 1, 1}, {1, 1, 1, 1}, {1, 1, 1, 1}, {1, 1, 1, 1}, {1, 1, 1, 1}, {1, 1, 1, 1}};
    final int[][] positions = {{}, {1, 2, 3, 4}, {1, 2, 3, 4}, {1, 2, 3, 4}, {1, 2, 3, 4}, {1, 2, 3, 4}, {1, 2, 3, 4}, {1, 2, 3, 4}, {1, 2, 3, 4}, {1, 2, 3, 4}, {1, 2, 3, 4}, {1, 2, 3, 4}, {1, 2, 3, 4}, {1, 2, 3, 4}, {1, 2, 3, 4}, {1, 2, 3, 4}};
    final PacBio.HitMap hitMap = mapHits(contigs, positions);
    assertEquals(1, hitMap.size());
    final String actual = hitMap.get(1L).toString();
    assertTrue(actual, actual.startsWith("[HitCollection[Hit{C1, R1}, Hit{C2, R2}, Hit{C3, R3}, Hit{C4, R4}], HitCollection[Hit{C2, R1}, Hit{C3, R2}"));

  }

  private PacBio.HitMap mapHits(long[][] contigs, int[][] positions) {
    final GraphKmerAttribute graph = GraphMapCliTest.makeGraph(4, new String[]{}, new long[][]{});
    final List<List<ContigPosition>> lists = new ArrayList<>();
    for (int i = 0; i < contigs.length; ++i) {
      final List<ContigPosition> posList = new ArrayList<>();
      for (int j = 0; j < contigs[i].length; ++j) {
        posList.add(new ContigPosition(contigs[i][j], positions[i][j], graph));
      }
      lists.add(posList);
    }
    return PacBio.joinHits(lists, 2);
  }
  static final String LONG_CONTIG =       "AACTCCGGCAGTAGGACTCGAACCTACGACATCATGATTAACAGTCATGCGCTACTACCAACTGAGCTATGCCGGAATAATCGCGTGGCGACGTCCTACTCTCACAAAGG";
  static final String READ =         "AAAATAACTC-GGCAGTAGGACTCGAACCTACGACATCATGATTAACAGTCATGGGCTAC-ACCAACTGAGCTATG-CGGAATAATCGCGTGGCGACGTCCTACTCTCAC-AAGGACATATA";
  public void testAlignHits() {
    final GraphKmerAttribute graph = GraphMapCliTest.makeGraph(5, new String[]{LONG_CONTIG}, new long[][] {});
    final PacBio.HitCollection hc = new PacBio.HitCollection(1);
    hc.add(new PacBio.HitPosition(30, 34));
    final byte[] read = DnaUtils.encodeStringWithHyphen(READ);
    final PartialAlignment partialAlignment = PacBio.alignHitRightFixed(hc, graph, read);
    assertEquals(1, partialAlignment.getContig());
    assertEquals(0, partialAlignment.getContigStart());
    assertEquals(30, partialAlignment.getContigEnd());
    assertEquals(5, partialAlignment.getReadStart());
    assertEquals(34, partialAlignment.getReadEnd());
    assertEquals(2, partialAlignment.getAlignmentScore());
    final PartialAlignment partialAlignmentForward = PacBio.alignHit(hc, graph, read);
    assertEquals(1, partialAlignmentForward.getContig());
    assertEquals(16, partialAlignmentForward.getContigStart());
    assertEquals(graph.contigLength(1), partialAlignmentForward.getContigEnd());
    assertEquals(20, partialAlignmentForward.getReadStart());
    assertEquals(read.length - 7, partialAlignmentForward.getReadEnd());
    assertEquals(9, partialAlignmentForward.getAlignmentScore());
  }
  static final String SHORTER_READ =         "CTC-GGCAGTAGGACTCGAACCTACGACATCATGATTAACAGTCATGGGCTAC-ACCAACTGAGCTATG-CGGAATAATCGCGTGGCGACGTCCTACTCTCACA";
  public void testAlignHitsShortRead() {
    final GraphKmerAttribute graph = GraphMapCliTest.makeGraph(5, new String[]{LONG_CONTIG}, new long[][] {});
    final PacBio.HitCollection hc = new PacBio.HitCollection(1);
    hc.add(new PacBio.HitPosition(30, 27));
    final byte[] read = DnaUtils.encodeStringWithHyphen(SHORTER_READ);
    final PartialAlignment partialAlignment = PacBio.alignHitRightFixed(hc, graph, read);
    assertEquals(1, partialAlignment.getContig());
    assertEquals(2, partialAlignment.getContigStart());
    assertEquals(30, partialAlignment.getContigEnd());
    assertEquals(0, partialAlignment.getReadStart());
    assertEquals(27, partialAlignment.getReadEnd());
    assertEquals(2, partialAlignment.getAlignmentScore());
    final PartialAlignment partialAlignmentForward = PacBio.alignHit(hc, graph, read);

    assertEquals(1, partialAlignmentForward.getContig());
    assertEquals(16, partialAlignmentForward.getContigStart());
    assertEquals(graph.contigLength(1) - 4, partialAlignmentForward.getContigEnd());
    assertEquals(13, partialAlignmentForward.getReadStart());
    assertEquals(read.length, partialAlignmentForward.getReadEnd());
    assertEquals(7, partialAlignmentForward.getAlignmentScore());
  }
  static final String MANY_INSERTS =     "AAACTC-GG-AGT-GGA-TCG-ACCTACGACATCATGATTAACAGTCATGGGCTAC-ACCAACTGA-CTATG-CGGAAT-ATCGCGT-GCGACGT-CTACT-TCAC-AAGGATAA";
  public void testAlignHitsHeapsOfInserts() {
    final GraphKmerAttribute graph = GraphMapCliTest.makeGraph(5, new String[]{LONG_CONTIG}, new long[][] {});
    final PacBio.HitCollection hc = new PacBio.HitCollection(1);
    hc.add(new PacBio.HitPosition(38, 34));
    final byte[] read = DnaUtils.encodeStringWithHyphen(MANY_INSERTS);
    final PartialAlignment partialAlignment = PacBio.alignHitRightFixed(hc, graph, read);
    assertEquals(1, partialAlignment.getContig());
    assertEquals(0, partialAlignment.getContigStart());
    assertEquals(38, partialAlignment.getContigEnd());
    assertEquals(1, partialAlignment.getReadStart());
    assertEquals(34, partialAlignment.getReadEnd());
    assertEquals(10, partialAlignment.getAlignmentScore());
    final PartialAlignment partialAlignmentForward = PacBio.alignHit(hc, graph, read);

    assertEquals(1, partialAlignmentForward.getContig());
    assertEquals(24, partialAlignmentForward.getContigStart());
    assertEquals(graph.contigLength(1), partialAlignmentForward.getContigEnd());
    assertEquals(20, partialAlignmentForward.getReadStart());
    assertEquals(read.length - 4, partialAlignmentForward.getReadEnd());
    assertEquals(19, partialAlignmentForward.getAlignmentScore());
  }

  static final String MANY_DELETES_READ =       "---TCCGGGCAGTAAGGAACTCGAAACCTAACGACATCATGATTAACAGTCATGCCGCTACCTACCAACTGAGGCTATGCCCGGAATAAATCGCCGTGGCGACGTCCTACTCTCAC";
  static final String MANY_DELETES_CONT =       "AACTCCGG-CAGT-AGGA-CTCG-AACCT-ACGACATCATGATTAACAGTCATG-CGCTA-CTACCAACTGAG-CTATGCC-GGAATAA-TCGC-GTGGCGACGTCCTACTCTCACAAAGG".replaceAll("-", "");
  public void testAlignHitsHeapsOfDeletes() {
    final GraphKmerAttribute graph = GraphMapCliTest.makeGraph(5, new String[]{MANY_DELETES_CONT}, new long[][] {});
    final PacBio.HitCollection hc = new PacBio.HitCollection(1);
    hc.add(new PacBio.HitPosition(43, 45));
    final byte[] read = DnaUtils.encodeStringWithHyphen(MANY_DELETES_READ);
    final PartialAlignment partialAlignment = PacBio.alignHitRightFixed(hc, graph, read);
    assertEquals(1, partialAlignment.getContig());
    assertEquals(3, partialAlignment.getContigStart());
    assertEquals(43, partialAlignment.getContigEnd());
    assertEquals(0, partialAlignment.getReadStart());
    assertEquals(45, partialAlignment.getReadEnd());
    assertEquals(10, partialAlignment.getAlignmentScore());
    final PartialAlignment partialAlignmentForward = PacBio.alignHit(hc, graph, read);

    assertEquals(1, partialAlignmentForward.getContig());
    assertEquals(29, partialAlignmentForward.getContigStart());
    assertEquals(graph.contigLength(1) - 5, partialAlignmentForward.getContigEnd());
    assertEquals(31, partialAlignmentForward.getReadStart());
    assertEquals(read.length, partialAlignmentForward.getReadEnd());
    assertEquals(12, partialAlignmentForward.getAlignmentScore());
  }




  static final String LONG_CONTIG_INSERT =       "GACTCC-GGCNGTAGGACTCGAACCTACGACATCATGATTAACAGTCATGCGCTAC-TACCAACTGAGCTATGC-GGAATAATCGCGTGGCGACGTCCTACTCTCACAA-GG".replaceAll("-", "");
  static final String READ_INSERT =         "AAAATGACTCCCGGCAGTAGGACTCGAACCTACGACATCATGATTAACAGTCATGGGCTACCTACCAACTGAGCTATGCCGGAATAATCGCGTGGCGACGTCCTACTCTCACAAAGGACATATA";
  public void testAlignHitsInserts() {
    final GraphKmerAttribute graph = GraphMapCliTest.makeGraph(5, new String[] {LONG_CONTIG_INSERT}, new long[][] {});
    final PacBio.HitCollection hc = new PacBio.HitCollection(1);
    hc.add(new PacBio.HitPosition(29, 35));
    final byte[] read = DnaUtils.encodeStringWithHyphen(READ_INSERT);
    final PartialAlignment partialAlignment = PacBio.alignHitRightFixed(hc, graph, read);
    assertEquals(1, partialAlignment.getContig());
    assertEquals(0, partialAlignment.getContigStart());
    assertEquals(29, partialAlignment.getContigEnd());
    assertEquals(5, partialAlignment.getReadStart());
    assertEquals(35, partialAlignment.getReadEnd());
    assertEquals(2, partialAlignment.getAlignmentScore());
    final PartialAlignment partialAlignmentForward = PacBio.alignHit(hc, graph, read);
    assertEquals(1, partialAlignmentForward.getContig());
    assertEquals(15, partialAlignmentForward.getContigStart());
    assertEquals(graph.contigLength(1), partialAlignmentForward.getContigEnd());
    assertEquals(21, partialAlignmentForward.getReadStart());
    assertEquals(read.length - 7, partialAlignmentForward.getReadEnd());
    assertEquals(9, partialAlignmentForward.getAlignmentScore());
  }
  static final String SHORTER_READ_INSERT =         "CTCCCGGCAGTAGGACTCGAACCTACGACATCATGATTAACAGTCATGGGCTACCTACCAACTGAGCTATGCCGGAATAATCGCGTGGCGACGTCCTACTCTCACA";
  public void testAlignHitsShortReadInserts() {
    final GraphKmerAttribute graph = GraphMapCliTest.makeGraph(5, new String[]{LONG_CONTIG_INSERT}, new long[][] {});
    final PacBio.HitCollection hc = new PacBio.HitCollection(1);
    hc.add(new PacBio.HitPosition(27, 26));
    final byte[] read = DnaUtils.encodeStringWithHyphen(SHORTER_READ_INSERT);
    final PartialAlignment partialAlignment = PacBio.alignHitRightFixed(hc, graph, read);
    assertEquals(1, partialAlignment.getContig());
    assertEquals(2, partialAlignment.getContigStart());
    assertEquals(27, partialAlignment.getContigEnd());
    assertEquals(0, partialAlignment.getReadStart());
    assertEquals(26, partialAlignment.getReadEnd());
    assertEquals(2, partialAlignment.getAlignmentScore());
    final PartialAlignment partialAlignmentForward = PacBio.alignHit(hc, graph, read);

    assertEquals(1, partialAlignmentForward.getContig());
    assertEquals(13, partialAlignmentForward.getContigStart());
    assertEquals(graph.contigLength(1) - 3, partialAlignmentForward.getContigEnd());
    assertEquals(12, partialAlignmentForward.getReadStart());
    assertEquals(read.length, partialAlignmentForward.getReadEnd());
    assertEquals(7, partialAlignmentForward.getAlignmentScore());
  }

  public void testJoinAlignments() {
    final GraphKmerAttribute graph = GraphMapCliTest.makeGraph(5, new String[]{"ACGTGTGTG", "ACCCACACT", "ACCTTAGATAGT", "TTGTTGGA", "GGGGGAGATTTTGAG"}, new long[][] {{1, 2}, {1, 3}, {4, 3}, {3, 5}});
    final List<PartialAlignment> partialAlignments = new ArrayList<>();
    final int part1 = graph.contigLength(1) - 1;
    final int part2Start = part1 - 5;
    final int part2End = part2Start + graph.contigLength(2);
    final int alt2Start = part1 - 5;
    final int alt2End = alt2Start + graph.contigLength(3);
    final int part3Start = alt2End - 5;
    final int part3End = part3Start + graph.contigLength(5) - 2;



    partialAlignments.add(new PartialAlignment(1, 0, part1, 1, 1, graph.contigLength(1)));
    partialAlignments.add(new PartialAlignment(1, part2Start, part2End, 2, 0, graph.contigLength(2)));
    partialAlignments.add(new PartialAlignment(1, alt2Start, alt2End, 3, 0, graph.contigLength(3)));
    partialAlignments.add(new PartialAlignment(1, part3Start, part3End, 5, 0, graph.contigLength(5) - 2));
    final Map<Long, List<PacBioPath>> longListMap = PacBio.joinAlignments(partialAlignments, graph);
    assertEquals(1L, longListMap.get(2L).get(0).mPrevious.mAlignment.getContig());
    assertEquals(2L, longListMap.get(2L).get(0).mAlignment.getContig());

    assertEquals(1L, longListMap.get(3L).get(0).mPrevious.mAlignment.getContig());
    assertEquals(1L, longListMap.get(5L).get(0).mPrevious.mPrevious.mAlignment.getContig());
  }

  public void testJoinAlignmentsBetter() {

    final GraphKmerAttribute graph = GraphMapCliTest.makeGraph(5, new String[]{"ACGTGTGTG", "ACCCACACT", "ACCTTAGATAGT", "TTGTTGGA", "GGGGGAGATTTTGAG"}, new long[][] {{1, 2}, {1, 3}, {4, 3}, {3, 5}});
    final List<PartialAlignment> partialAlignments = new ArrayList<>();
    final int part1 = graph.contigLength(1) - 1;
    final int part2Start = part1 - 5;
    final int part2End = part2Start + graph.contigLength(2);
    final int alt2Start = part1 - 5;
    final int alt2End = alt2Start + graph.contigLength(3);
    final int part3Start = alt2End - 5;
    final int part3End = part3Start + graph.contigLength(5) - 2;



    partialAlignments.add(new PartialAlignment(3, 0, part1, 1, 1, graph.contigLength(1)));
    partialAlignments.add(new PartialAlignment(1, 0, part1, 1, 1, graph.contigLength(1)));
    partialAlignments.add(new PartialAlignment(1, part2Start, part2End, 2, 0, graph.contigLength(2)));
    partialAlignments.add(new PartialAlignment(1, alt2Start, alt2End, 3, 0, graph.contigLength(3)));
    partialAlignments.add(new PartialAlignment(1, part3Start, part3End, 5, 0, graph.contigLength(5) - 2));
    final Map<Long, List<PacBioPath>> longListMap = PacBio.joinAlignments(partialAlignments, graph);
    assertEquals(1L, longListMap.get(2L).get(0).mPrevious.mAlignment.getContig());
    assertEquals(2L, longListMap.get(2L).get(0).mAlignment.getContig());

    assertEquals(1L, longListMap.get(3L).get(0).mPrevious.mAlignment.getContig());
    assertEquals(1L, longListMap.get(5L).get(0).mPrevious.mPrevious.mAlignment.getContig());
    assertEquals(3, longListMap.get(5L).get(0).mScore);
  }
  public void testJoinAlignmentsReplace() {
    final GraphKmerAttribute graph = GraphMapCliTest.makeGraph(5, new String[]{"ACGTGTGTG", "ACCCACACT", "ACCTTAGAT", "TTGTTGGA"}, new long[][] {{1, 2}, {1, 3}, {3, 4}, {2, 4}});
    final List<PartialAlignment> partialAlignments = new ArrayList<>();
    final int part1 = graph.contigLength(1) - 1;
    final int part2Start = part1 - 5;
    final int part2End = part2Start + graph.contigLength(2);
    final int part3Start = part2End - 5;
    final int part3End = part3Start + graph.contigLength(4) - 2;
    partialAlignments.add(new PartialAlignment(3, 0, part1, 1, 1, graph.contigLength(1)));
    partialAlignments.add(new PartialAlignment(1, 0, part1, 1, 1, graph.contigLength(1)));
    partialAlignments.add(new PartialAlignment(2, part2Start, part2End - 1, 2, 0, graph.contigLength(2)));
    partialAlignments.add(new PartialAlignment(2, part2Start, part2End, 2, 0, graph.contigLength(2)));
    partialAlignments.add(new PartialAlignment(1, part2Start, part2End, 3, 0, graph.contigLength(3)));
    partialAlignments.add(new PartialAlignment(1, part3Start, part3End, 4, 0, graph.contigLength(4) - 2));

    final Map<Long, List<PacBioPath>> longListMap = PacBio.joinAlignments(partialAlignments, graph);
    assertEquals(1L, longListMap.get(2L).get(0).mPrevious.mAlignment.getContig());
    assertEquals(2L, longListMap.get(2L).get(0).mAlignment.getContig());

    assertEquals(1L, longListMap.get(3L).get(0).mPrevious.mAlignment.getContig());
    assertEquals(1L, longListMap.get(4L).get(0).mPrevious.mPrevious.mAlignment.getContig());
    assertEquals(1, longListMap.get(4L).size());
    assertEquals(3L, longListMap.get(4L).get(0).mPrevious.mAlignment.getContig());
    assertEquals(3, longListMap.get(4L).get(0).score());
  }
  public void testIsInternal() {
    final Graph g = GraphMapCliTest.makeGraph(4, new String[] {"AAACCGGAG", "AGATAGATGGTA"}, new long[][] {});
    final List<List<ContigPosition>> hits = new ArrayList<>();
    addEmptyHits(hits, 3);
    hits.add(Arrays.asList(new ContigPosition(2, 3, g)));
    hits.add(Arrays.asList(new ContigPosition(2, 4, g), new ContigPosition(1, 1, g)));
    hits.add(Arrays.asList(new ContigPosition(2, 5, g)));
    addEmptyHits(hits, 3);
    hits.add(Arrays.asList(new ContigPosition(2, 9, g)));
    addEmptyHits(hits, 2);
    assertTrue(PacBio.isInternal(g, hits));
    // Make read longer than contig
    addEmptyHits(hits, 1);
    assertFalse(PacBio.isInternal(g, hits));
  }
  public void testIsInternalOverlapStart() {
    final Graph g = GraphMapCliTest.makeGraph(4, new String[] {"AAACCGGAG", "AGATAGATGGTA"}, new long[][] {});
    final List<List<ContigPosition>> hits = new ArrayList<>();
    addEmptyHits(hits, 4);
    hits.add(Arrays.asList(new ContigPosition(2, 3, g)));
    assertFalse(PacBio.isInternal(g, hits));
  }
  public void testIsInternalOverlapEnd() {
    final Graph g = GraphMapCliTest.makeGraph(4, new String[] {"AAACCGGAG", "AGATAGATGGTA"}, new long[][] {});
    final List<List<ContigPosition>> hits = new ArrayList<>();
    addEmptyHits(hits, 4);
    hits.add(Arrays.asList(new ContigPosition(2, 8, g)));
    addEmptyHits(hits, 3);
    assertTrue(PacBio.isInternal(g, hits));
    // Make read overlap end
    addEmptyHits(hits, 1);
    assertFalse(PacBio.isInternal(g, hits));
  }

  public void testNoHits() {
    final Graph g = GraphMapCliTest.makeGraph(4, new String[] {"AAACCGGAG", "AGATAGATGGTA"}, new long[][] {});
    final List<List<ContigPosition>> hits = new ArrayList<>();
    addEmptyHits(hits, 4);
    assertFalse(PacBio.isInternal(g, hits));
  }

  private void addEmptyHits(List<List<ContigPosition>> hits, int emptyCount) {
    for (int i = 0; i < emptyCount; ++i) {
      hits.add(Collections.<ContigPosition>emptyList());
    }
  }

  public void testBestPath() {
    final PacBioPath first = new PacBioPath(null, new PartialAlignment(2, 0, 20, 1, 11, 21));
    final PacBioPath second = new PacBioPath(first, new PartialAlignment(1, 15, 25, 2, 1, 11));
    final PacBioPath alternate = new PacBioPath(first, new PartialAlignment(3, 15, 25, 3, 1, 11));
    final PacBioPath sameScore = new PacBioPath(first, new PartialAlignment(1, 15, 25, 4, 1, 11));
    final Map<Long, List<PacBioPath>> paths = new HashMap<>();
    paths.put(1L, Arrays.asList(first));
    paths.put(2L, Arrays.asList(second));
    paths.put(3L, Arrays.asList(alternate));
    assertEquals(null, PacBio.uniqueBest(paths, 26));
    assertEquals(second, PacBio.uniqueBest(paths, 25));
    paths.put(4L, Arrays.asList(sameScore));
    assertEquals(null, PacBio.uniqueBest(paths, 25));
  }

  public void testEndToEnd() throws IOException {
    final ReadPairSource source = new ReadPairSource(ReaderTestUtils.getReaderDnaMemory(ReaderTestUtils.fasta(
        "GAGATATATGATAGATAGACCCCACATAGATACAGGAGGGATTTAGAGGATAGGAGATGGATTGTGGAGCCGTGCCCGCC")));
    final String[] contigs = {
        "ACCGAGATATATGATAGATAGACCCAC"
                    , "GATAGATAGACCCACATAGATACAGGAGGGATTTGAGGATGGAGATGG"
                                                          , "TGAGGATGGAGATGGATTGTGGAGCCGTGCCCGCCCCCAAG"
                                                          , "TGAGGATGGAGATGGATTTTATTACCCCCACACTGGGGGAT"};
    final long[][] paths = {{1, 2}, {2, 3}, {2, 4}};
    final MutableGraph graph = GraphMapCliTest.makeGraph(15, contigs, paths);
    new PacBio(PacBioParams.builder().create(), TestUtils.getNullOutputStream()).mapPacBio(source, graph);
    assertEquals(4, graph.numberPaths());
    assertEquals(Arrays.asList(1L, 2L, 3L), PathArray.toList(graph.path(-4)));
  }

  public void testPrint() {
    // Test insert/delete
    final MemoryPrintStream mps = new MemoryPrintStream();
    PacBio.printAlignment("foo", new byte[] {1, 2, 3, 4, 1, 3, 4}, new byte[] {1, 2, 4, 1, 2, 3, 4}, 1, 1, "=D==I=", mps.printStream());
    final String expected = "foo" + LS
        + "C-TACG" + LS
        + "=D==I=" + LS
        + "CGTA-G" + LS
        + LS;
    assertEquals(expected, mps.toString());


    // Test line wrapping
    final byte[] contig = new byte[160];
    final byte[] read = new byte[160];
    final StringBuilder builder = new StringBuilder();
    for (int i = 0; i < 160; ++i) {
      contig[i] = 1;
      read[i] = 2;
      builder.append("X");
    }
    mps.reset();
    PacBio.printAlignment("foo", read, contig, 0, 0, builder.toString(), mps.printStream());
    final StringBuilder readString = new StringBuilder();
    final StringBuilder contigString = new StringBuilder();
    final StringBuilder match = new StringBuilder();
    for (int i = 0; i < 120; ++i) {
      contigString.append("A");
      readString.append("C");
      match.append("X");
    }
    assertTrue(mps.toString().contains(LS + contigString + LS + match + LS + readString + LS));
  }
}
