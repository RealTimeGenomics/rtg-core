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
package com.rtg.variant.bayes.snp;

import static com.rtg.util.StringUtils.LS;

import java.io.ByteArrayOutputStream;
import java.io.Closeable;
import java.io.PrintStream;

import com.rtg.mode.DnaUtils;
import com.rtg.sam.BadSuperCigarException;
import com.rtg.sam.SamUtils;
import com.rtg.util.MathUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.machine.MachineType;
import com.rtg.variant.MachineErrorParams;
import com.rtg.variant.VariantAlignmentRecord;
import com.rtg.variant.VariantParams;
import com.rtg.variant.VariantParamsBuilder;
import com.rtg.variant.bayes.EvidenceInterface;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import junit.framework.TestCase;

/**
 */
public class CigarParserModelTest extends TestCase {

  @Override
  protected void setUp() {
    Diagnostic.setLogStream();
  }

  private static class MockModelMatcher implements Closeable, MatcherInterface {
    private final ByteArrayOutputStream mBOS;
    private final PrintStream mPR;

    private final ByteArrayOutputStream mInsertBOS;
    private final PrintStream mInsertPR;
    private int mCount = 0;

    MockModelMatcher() {
      mBOS = new ByteArrayOutputStream();
      mPR = new PrintStream(mBOS);
      mInsertBOS = new ByteArrayOutputStream();
      mInsertPR = new PrintStream(mInsertBOS);
    }

    static int probToScore(final double p) {
      return (int) MathUtils.round(-10.0 * Math.log10(p));
    }

    @Override
    public void match(final int ref, final EvidenceInterface ev) {
      if (ev != null) {
        if (ev instanceof  EvidenceIndel) {
          mInsertPR.print("ref=" + ref);
          mInsertPR.println(" n");
        } else {
          mPR.print("ref=" + ref);
          mPR.print(" nt=" + ev.read() + " q=" + probToScore(ev.error()));
          mPR.println(" r=" + probToScore(ev.mapError()));
          mCount++;
        }
      }
    }

    int count() {
      return mCount;
    }

    @Override
    public void match(int refPosition, int readBasesLeft, int readBasesRight, int readNt, int mapQ, int phred, int stateIndex) {
      mPR.println("ref=" + refPosition + " nt=" + readNt + " q="
          + phred + " r=" + mapQ);
      mCount++;
    }

    @Override
    public void unmapped(int refPosition) {
    }

    @Override
    public int getStateIndex(boolean isForward, boolean isReadPaired, boolean isMated) {
      return 0;
    }

    @Override
    public void close() {
      mPR.close();
      mInsertPR.close();
    }

    @Override
    public String toString() {
      return mBOS.toString();
    }

    /**
     * Equivalent of toString but displays the inserts
     * @return string representation of the insert matches
     */
    public String insertString() {
      return mInsertBOS.toString();
    }
  }

  private static VariantAlignmentRecord toVariantAlignmentRecord(byte[] read, String qualities, String cigar, int zeroBasedStart, boolean forward) {
    final SAMRecord sam = new SAMRecord(new SAMFileHeader());
    sam.setCigarString(cigar);
    sam.setReadBases(read);
    sam.setAlignmentStart(zeroBasedStart + 1);
    sam.setBaseQualityString(qualities == null ? "*" : qualities);
    if (forward) {
      sam.setFlags(67); //first in pair, not reverse
    } else {
      sam.setFlags(131);  //second in pair, not reverse.
    }

    return new VariantAlignmentRecord(sam);
  }

  public void testCigarToMatcherOffTemplate() {
    final MachineErrorParams me = MachineErrorParams.builder().create();

    try (MockModelMatcher matcher = new MockModelMatcher()) {
      final String seq = "GAGAAGCTGGAGGCGCCGGCAACCAGCTGGGGATA";
      final byte[] templateBytes = DnaUtils.encodeString("AAGAGAAGCTGGTTTTTTAGGCGCCGGCAACCAGCTGGAT");
      final VariantAlignmentRecord var = toVariantAlignmentRecord(seq.getBytes(), null, "10M6N20M2I3M", 2, true);
      final CigarParserModel cpm = new CigarParserModel(matcher, null, 0, templateBytes.length, new VariantParamsBuilder().create()) {
        @Override
        int getReadScore(VariantAlignmentRecord var) {
          return 63;
        }
      };
      cpm.toMatcher(var, me.machineType(), 20, templateBytes);
      fail();
    } catch (final BadSuperCigarException e) {
      // Expected
    }
    // assertEquals(CIGAR_EXP1, matcher.toString());
  }
  private static final String CIGAR_EXP1 = ""
      + "ref=2 nt=1 q=63 r=63" + LS
      + "ref=3 nt=2 q=0 r=63" + LS
      + "ref=4 nt=2 q=63 r=63" + LS
      + "ref=5 nt=3 q=0 r=63" + LS
      ;

  private static void check(String seq, byte[] templateBytes, int regionStart, int regionEnd, int templateStart, String cigar, String qual, String exp, String expInsert) throws Exception {
    final MockModelMatcher matcher = new MockModelMatcher();
    final VariantParams vp = new VariantParamsBuilder().ignoreQualityScores(false).create();
    try {
      final CigarParserModel cpm = new CigarParserModel(matcher, matcher, regionStart, regionEnd, vp) {
        @Override
        int getReadScore(VariantAlignmentRecord var) {
          return 63;
        }
      };
      final VariantAlignmentRecord var = toVariantAlignmentRecord(seq.getBytes(), qual, cigar, templateStart, true);
      cpm.toMatcher(var, null, 20, templateBytes);
    } finally {
      matcher.close();
    }
    assertEquals(exp, matcher.toString());
    if (expInsert != null) {
      assertEquals(expInsert, matcher.insertString());
    }
  }

  private static final byte[] TEMPLATEBYTES = {0, 1, 1, 2, 2, 3, 4};
  private static final byte[] TEMPLATEBYTES2 = {0, 0, 2, 3, 4, 1, 2};

  // matches
  public void testCigarToMatcher1() throws Exception {

    final String seq = "ACCG";

    check(seq, TEMPLATEBYTES, 0, TEMPLATEBYTES.length, 2, "4M", "`!`!", CIGAR_EXP1, null);
    check(seq, TEMPLATEBYTES, 0, TEMPLATEBYTES.length, 2, "4=", "`!`!", CIGAR_EXP1, null);
    check(seq, TEMPLATEBYTES2, 0, TEMPLATEBYTES2.length, 2, "4X", "`!`!", CIGAR_EXP1, null);

  }

  private static final String CIGAR_EXP2 = ""
      + "ref=1 nt=1 q=63 r=63" + LS
      + "ref=2 nt=1 q=0 r=63" + LS
      + "ref=5 nt=3 q=63 r=63" + LS
      + "ref=6 nt=4 q=0 r=63" + LS
      ;

  // deletion
  public void testCigarToMatcher2() throws Exception {
    final String seq = "AAGT";

    check(seq, TEMPLATEBYTES, 0, TEMPLATEBYTES.length, 1, "2M2D2M", "`!`!", CIGAR_EXP2, null);
    check(seq, TEMPLATEBYTES, 0, TEMPLATEBYTES.length, 1, "2=2D2=", "`!`!", CIGAR_EXP2, null);
    check(seq, TEMPLATEBYTES2, 0, TEMPLATEBYTES2.length, 1, "2X2D2X", "`!`!", CIGAR_EXP2, null);

  }

  public void testCigarToMatcherBad1() {
    MockModelMatcher matcher = new MockModelMatcher();
    final VariantParams vp = new VariantParamsBuilder().ignoreQualityScores(false).create();

    final String seq = "ACRG";
    final byte[] templateBytes = {0, 0, 1, 2, 2, 3, 4};
    try {
      final CigarParserModel cpm = new CigarParserModel(matcher, matcher, 0, templateBytes.length, vp) {
        @Override
        int getReadScore(VariantAlignmentRecord var) {
          return 63;
        }
      };
      final VariantAlignmentRecord var = toVariantAlignmentRecord(seq.getBytes(), "`!`!", "2M2D1N1M", 1, true);
      cpm.toMatcher(var, null, 20, templateBytes);
      fail();
    } catch (final BadSuperCigarException e) {
      // Expected
    } finally {
      matcher.close();
    }

    matcher = new MockModelMatcher();
    try {
      final CigarParserModel cpm = new CigarParserModel(matcher, matcher, 0, templateBytes.length, vp) {
        @Override
        int getReadScore(VariantAlignmentRecord var) {
          return 63;
        }
      };
      final VariantAlignmentRecord var = toVariantAlignmentRecord(seq.getBytes(), "`!`!", "2=2D1N1=", 1, true);
      cpm.toMatcher(var, null, 20, templateBytes);
      fail();
    } catch (final BadSuperCigarException e) {
      // Expected
    } finally {
      matcher.close();
    }

    matcher = new MockModelMatcher();
    try {
      final CigarParserModel cpm = new CigarParserModel(matcher, matcher, 0, templateBytes.length, vp) {
        @Override
        int getReadScore(VariantAlignmentRecord var) {
          return 63;
        }
      };
      final VariantAlignmentRecord var = toVariantAlignmentRecord(seq.getBytes(), "`!`!", "2X2D1N1X", 1, true);
      cpm.toMatcher(var, null, 20, templateBytes);
      fail();
    } catch (final BadSuperCigarException e) {
      // Expected
    } finally {
      matcher.close();
    }
  }

  private static final String ERR2 = "Malformed CIGAR string: 2MDN7M";

  public void testCigarToMatcherBad2() throws Exception {
    final VariantParams vp = new VariantParamsBuilder().ignoreQualityScores(false).create();
    try (MockModelMatcher matcher = new MockModelMatcher()) {
      final String seq = "ACTG";
      final byte[] templateBytes = {0, 0, 1, 2, 2, 3, 4};
      final CigarParserModel cpm = new CigarParserModel(matcher, matcher, 0, templateBytes.length, vp) {
        @Override
        int getReadScore(VariantAlignmentRecord var) {
          return 63;
        }
      };
      final VariantAlignmentRecord var = toVariantAlignmentRecord(seq.getBytes(), "`!`!", "2MDN7M", 1, true);
      cpm.toMatcher(var, null, 20, templateBytes);
      fail();
    } catch (final RuntimeException e) {
      assertEquals(ERR2, e.getMessage());
    }
  }

  private static final String CIGAR_EXP3 = ""
      + "ref=1 nt=4 q=63 r=63" + LS
      + "ref=2 nt=3 q=0 r=63" + LS;
  private static final String INSERT_EXP3 = ""
      + "ref=1 n" + LS
      ;

  public void testInsertion() throws Exception {
    final String seq = "ACTG";
    check(seq, TEMPLATEBYTES, 0, TEMPLATEBYTES.length, 1, "2I2M", "`!`!", CIGAR_EXP3, INSERT_EXP3);
  }

  private static final String CIGAR_EXP3R = ""
      + "ref=2 nt=3 q=0 r=63" + LS;

  public void testInsertionR() throws Exception {
    final String seq = "ACTG";
    check(seq, TEMPLATEBYTES, 2, TEMPLATEBYTES.length, 1, "2I2M", "`!`!", CIGAR_EXP3R, "");
  }

  private static final String CIGAR_EXP_AT_END = ""
      + "ref=1 nt=1 q=63 r=63" + LS
      + "ref=2 nt=2 q=0 r=63" + LS
      ;
  private static final String INSERT_EXP_AT_END = ""
      + "ref=3 n" + LS
      ;

  public void testInsertionAtEnd() throws Exception {
    final String seq = "ACTG";
    check(seq, TEMPLATEBYTES, 0, TEMPLATEBYTES.length, 1, "2M2I", "`!`!", CIGAR_EXP_AT_END, INSERT_EXP_AT_END);
  }

  private static final String CIGAR_EXP_AT_ENDR = ""
      + "ref=1 nt=1 q=63 r=63" + LS
      ;

  public void testInsertionAtEndR() throws Exception {
    final String seq = "ACTG";
    check(seq, TEMPLATEBYTES, 1, 2, 1, "2M2I", "`!`!", CIGAR_EXP_AT_ENDR, "");
  }

  private static final String CIGAR_EXP4 = ""
      + "ref=2 nt=1 q=0 r=63" + LS
      + "ref=3 nt=2 q=63 r=63" + LS
      ;

  // insertion
  public void testCigarToMatcher4() throws Exception {
    final String seq = "TACT";

    check(seq, TEMPLATEBYTES, 0, TEMPLATEBYTES.length, 1, "1I1D2M1I", "`!`!", CIGAR_EXP4, null);
    check(seq, TEMPLATEBYTES, 0, TEMPLATEBYTES.length, 1, "1I1D2=1I", "`!`!", CIGAR_EXP4, null);
    check(seq, TEMPLATEBYTES2, 0, TEMPLATEBYTES2.length, 1, "1I1D2X1I", "`!`!", CIGAR_EXP4, null);
  }

  private static final String CIGAR_NS = "" + "ref=1 nt=1 q=63 r=63" + LS + "ref=4 nt=2 q=63 r=63" + LS;

  // deletion
  public void testCigarToMatcherNs() throws Exception {
    final String seq = "AC";
    final byte[] templateBytes = {0, 1, 0, 0, 2};
    final byte[] templateBytes2 = {0, 2, 0, 0, 1};

    check(seq, templateBytes, 0, templateBytes.length, 1, "1M2N1M", "``", CIGAR_NS, null);
    check(seq, templateBytes, 0, templateBytes.length, 1, "1=2N1=", "``", CIGAR_NS, null);
    check(seq, templateBytes2, 0, templateBytes2.length, 1, "1X2N1X", "``", CIGAR_NS, null);

  }

  private static final String CIGAR_EXP5 = ""
      + "ref=3 nt=2 q=63 r=63" + LS
      + "ref=4 nt=2 q=0 r=63" + LS
      + "ref=5 nt=3 q=63 r=63" + LS
      + "ref=6 nt=4 q=0 r=63" + LS
      ;

  // deletion
  public void testCigarToMatcher5() throws Exception {
    final String seq = "CCGT";
    check(seq, TEMPLATEBYTES, 0, TEMPLATEBYTES.length, 1, "2D4M", "`!`!", CIGAR_EXP5, null);
    check(seq, TEMPLATEBYTES, 0, TEMPLATEBYTES.length, 1, "2D4=", "`!`!", CIGAR_EXP5, null);
    check(seq, TEMPLATEBYTES2, 0, TEMPLATEBYTES2.length, 1, "2D4X", "`!`!", CIGAR_EXP5, null);
  }

  private static final String CIGAR_EXP6 = ""
      + "ref=1 nt=1 q=63 r=63" + LS
      + "ref=2 nt=1 q=0 r=63" + LS
      + "ref=3 nt=2 q=63 r=63" + LS
      + "ref=4 nt=2 q=0 r=63" + LS
      ;

  private static final String CIGAR_EXP6_DEL = ""
      + "ref=5 n" + LS
      + "ref=6 n" + LS
      ;


  // deletion
  public void testCigarToMatcher6() throws Exception {
    final String seq = "AACC";

    check(seq, TEMPLATEBYTES, 0, TEMPLATEBYTES.length, 1, "4M2D", "`!`!", CIGAR_EXP6, CIGAR_EXP6_DEL);
    check(seq, TEMPLATEBYTES, 0, TEMPLATEBYTES.length, 1, "4=2D", "`!`!", CIGAR_EXP6, CIGAR_EXP6_DEL);
    check(seq, TEMPLATEBYTES2, 0, TEMPLATEBYTES2.length, 1, "4X2D", "`!`!", CIGAR_EXP6, CIGAR_EXP6_DEL);

  }

  private static final String CIGAR_EXP6R = ""
      + "ref=3 nt=2 q=63 r=63" + LS
      + "ref=4 nt=2 q=0 r=63" + LS
      ;

  private static final String CIGAR_EXP6R_DEL = ""
      + "ref=5 n" + LS
      ;

  // deletion
  public void testCigarToMatcher6R() throws Exception {
    final String seq = "AACC";
    check(seq, TEMPLATEBYTES, 3, 6, 1, "4M2D", "`!`!", CIGAR_EXP6R, CIGAR_EXP6R_DEL);
    check(seq, TEMPLATEBYTES, 3, 6, 1, "4=2D", "`!`!", CIGAR_EXP6R, CIGAR_EXP6R_DEL);
    check(seq, TEMPLATEBYTES2, 3, 6, 1, "4X2D", "`!`!", CIGAR_EXP6R, CIGAR_EXP6R_DEL);
  }

  private static final String CIGAR_EXP7 = ""
      + "ref=1 nt=1 q=63 r=63" + LS
      + "ref=2 nt=2 q=0 r=63" + LS
      + "ref=3 nt=4 q=63 r=63" + LS
      + "ref=4 nt=3 q=0 r=63" + LS
      + "ref=5 nt=1 q=63 r=63" + LS
      + "ref=6 nt=2 q=63 r=63" + LS
      + "ref=7 nt=4 q=63 r=63" + LS
      + "ref=8 nt=3 q=63 r=63" + LS
      + "ref=9 nt=1 q=63 r=63" + LS
      + "ref=10 nt=2 q=63 r=63" + LS
      + "ref=11 nt=3 q=63 r=63" + LS
      + "ref=12 nt=4 q=63 r=63" + LS;

  // match
  public void testCigarToMatcher7() throws Exception {
    final String seq = "ACTGACTGACGT";
    final byte[] templateBytes = DnaUtils.encodeString("NACTGACTGACGT");
    final byte[] templateBytes2 = DnaUtils.encodeString("NGTCAGTCAGTAC");
    check(seq, templateBytes, 0, templateBytes.length, 1, "12M", "`!`!````````", CIGAR_EXP7, null);
    check(seq, templateBytes, 0, templateBytes.length, 1, "12=", "`!`!````````", CIGAR_EXP7, null);
    check(seq, templateBytes2, 0, templateBytes2.length, 1, "12X", "`!`!````````", CIGAR_EXP7, null);
  }

  private static final String CIGAR_EXP7R = ""
      + "ref=3 nt=4 q=63 r=63" + LS
      + "ref=4 nt=3 q=0 r=63" + LS
      + "ref=5 nt=1 q=63 r=63" + LS
      + "ref=6 nt=2 q=63 r=63" + LS
      + "ref=7 nt=4 q=63 r=63" + LS
      ;

  // match - restricted region
  public void testCigarToMatcher7R() throws Exception {
    final String seq = "ACTGACTGACGT";
    final byte[] templateBytes = DnaUtils.encodeString("NACTGACTGACGT");
    final byte[] templateBytes2 = DnaUtils.encodeString("NGTCAGTCAGTAC");
    check(seq, templateBytes, 3, 8, 1, "12M", "`!`!````````", CIGAR_EXP7R, null);
    check(seq, templateBytes, 3, 8, 1, "12=", "`!`!````````", CIGAR_EXP7R, null);
    check(seq, templateBytes2, 3, 8, 1, "12X", "`!`!````````", CIGAR_EXP7R, null);
  }


  private static final String CIGAR_EXP8 = "" + "ref=0 nt=4 q=63 r=63" + LS + "ref=1 nt=3 q=0 r=63"
      + LS + "ref=2 nt=2 q=63 r=63" + LS + "ref=3 nt=3 q=0 r=63" + LS;

  // soft-clipping
  public void testCigarToMatcher8() throws Exception {
    final String seq = "ACTGCG";
    check(seq, TEMPLATEBYTES3, 0, TEMPLATEBYTES3.length, 0, "2S4M", "`!`!`!", CIGAR_EXP8, null);
    check(seq, TEMPLATEBYTES3, 0, TEMPLATEBYTES3.length, 0, "2S4=", "`!`!`!", CIGAR_EXP8, null);
    check(seq, TEMPLATEBYTES4, 0, TEMPLATEBYTES4.length, 0, "2S4X", "`!`!`!", CIGAR_EXP8, null);
  }

  private static final String CIGAR_EXP9 = "" + "ref=0 nt=4 q=20 r=63" + LS + "ref=1 nt=3 q=20 r=63"
      + LS + "ref=2 nt=2 q=20 r=63" + LS + "ref=3 nt=3 q=20 r=63" + LS;

  private static final byte[] TEMPLATEBYTES3 = {4, 3, 2, 3, 0, 0};
  private static final byte[] TEMPLATEBYTES4 = {1, 2, 3, 4, 0, 0};

  // default q-values
  public void testCigarToMatcher9() throws Exception {
    final String seq = "ACTGCG";
    check(seq, TEMPLATEBYTES3, 0, TEMPLATEBYTES3.length, 0, "2S4M", null, CIGAR_EXP9, null);
    check(seq, TEMPLATEBYTES3, 0, TEMPLATEBYTES3.length, 0, "2S4=", null, CIGAR_EXP9, null);
    check(seq, TEMPLATEBYTES4, 0, TEMPLATEBYTES4.length, 0, "2S4X", null, CIGAR_EXP9, null);
  }

  private static final String CIGAR_EXP10 = "" + "ref=3 nt=1 q=20 r=63" + LS + "ref=4 nt=2 q=20 r=63"
      + LS + "ref=5 nt=4 q=20 r=63" + LS + "ref=6 nt=3 q=20 r=63" + LS;

  // case when match extends beyond end of template
  public void testCigarToMatcher10Bad() {
    final VariantParams vp = new VariantParamsBuilder().ignoreQualityScores(false).create();
    MockModelMatcher matcher = new MockModelMatcher();
    final String seq = "ACTGCG";
    final byte[] templateBytes = {0, 0, 1, 2, 2, 3, 4};
    try {
      final CigarParserModel cpm = new CigarParserModel(matcher, matcher, 0, templateBytes.length, vp) {
        @Override
        int getReadScore(VariantAlignmentRecord var) {
          return 63;
        }
      };
      final VariantAlignmentRecord var = toVariantAlignmentRecord(seq.getBytes(), null, "6M", 3, true);
      cpm.toMatcher(var, null, 20, templateBytes);
      fail();
    } catch (final BadSuperCigarException e) {
      // Expected
    } finally {
      matcher.close();
    }
    assertEquals("", matcher.toString());

    matcher = new MockModelMatcher();
    try {
      final CigarParserModel cpm = new CigarParserModel(matcher, matcher, 0, templateBytes.length, vp) {
        @Override
        int getReadScore(VariantAlignmentRecord var) {
          return 63;
        }
      };
      final VariantAlignmentRecord var = toVariantAlignmentRecord(seq.getBytes(), null, "6=", 3, true);
      cpm.toMatcher(var, null, 20, templateBytes);
      fail();
    } catch (final BadSuperCigarException e) {
      // Expected
    } finally {
      matcher.close();
    }
    assertEquals("", matcher.toString());

    matcher = new MockModelMatcher();
    try {
      final CigarParserModel cpm = new CigarParserModel(matcher, matcher, 0, templateBytes.length, vp) {
        @Override
        int getReadScore(VariantAlignmentRecord var) {
          return 63;
        }
      };
      final VariantAlignmentRecord var = toVariantAlignmentRecord(seq.getBytes(), null, "6X", 3, true);
      cpm.toMatcher(var, null, 20, templateBytes);
      fail();
    } catch (final BadSuperCigarException e) {
      // Expected
    } finally {
      matcher.close();
    }
    assertEquals("", matcher.toString());
  }

  public void testCigarToMatcher10Good() throws Exception {
    final String seq = "ACTGCG";
    final byte[] templateBytes1 = {0, 0, 1, 1, 2, 4, 3};
    final byte[] templateBytes2 = {0, 0, 1, 2, 3, 1, 4};
    check(seq, templateBytes1, 0, templateBytes1.length, 3, "4M2S", null, CIGAR_EXP10, null);
    check(seq, templateBytes1, 0, templateBytes1.length, 3, "4=2S", null, CIGAR_EXP10, null);
    check(seq, templateBytes2, 0, templateBytes2.length, 3, "4X2S", null, CIGAR_EXP10, null);
  }

  public void testGetDNA() {
    assertEquals(0, CigarParserModel.getDNA('N'));
    assertEquals(1, CigarParserModel.getDNA('A'));
    assertEquals(2, CigarParserModel.getDNA('C'));
    assertEquals(3, CigarParserModel.getDNA('G'));
    assertEquals(4, CigarParserModel.getDNA('T'));
  }

  public void testBugFoundByStu() throws Exception {
    final MockModelMatcher matcher = new MockModelMatcher();
    final VariantParams vp = new VariantParamsBuilder().ignoreQualityScores(false).create();
    try {
      final String seq = "CATATTAGGGCCTAATTGTAATACA";
      final byte[] templateBytes = DnaUtils.encodeString("GCATATTAGGGCCTATAATTGTAATAC");
      final CigarParserModel cpm = new CigarParserModel(matcher, matcher, 0, templateBytes.length, vp) {
        @Override
        int getReadScore(VariantAlignmentRecord var) {
          return 63;
        }
      };
      final VariantAlignmentRecord var = toVariantAlignmentRecord(seq.getBytes(), null, "1S13M1D3M1D8M", 1, true);
      cpm.toMatcher(var, MachineType.ILLUMINA_PE, 20, templateBytes);
    } finally {
      matcher.close();
    }
    assertNotNull(matcher.toString());
  }

  public void testCgLeftNoPrune() throws Exception {
    final MockModelMatcher matcher = new MockModelMatcher();
    final VariantParams vp = new VariantParamsBuilder().ignoreQualityScores(false).create();
    try {
      final String seq = "CCCTTTTTTTTTTTTTTTTTTTTAAAAAAAAAA";
      final byte[] templateBytes = DnaUtils.encodeString("AACCCTTTTTTTTTTTTTTTTTTTTGGGGGAAAAAAAAAA");
      final CigarParserModel cpm = new CigarParserModel(matcher, matcher, 0, templateBytes.length, vp) {
        @Override
        int getReadScore(VariantAlignmentRecord var) {
          return 63;
        }
      };
      final VariantAlignmentRecord var = toVariantAlignmentRecord(seq.getBytes(), null, "23=5N10=", 2, true);
      cpm.toMatcher(var, MachineType.COMPLETE_GENOMICS, 20, templateBytes);
    } finally {
      matcher.close();
    }
    assertEquals(""
      + (MachineType.CG_TRIM ?  ""
      : "ref=2 nt=2 q=20 r=63" + LS
        + "ref=3 nt=2 q=20 r=63" + LS
        + "ref=4 nt=2 q=20 r=63" + LS)
        + "ref=5 nt=4 q=20 r=63" + LS
        + "ref=6 nt=4 q=20 r=63" + LS
        + "ref=7 nt=4 q=20 r=63" + LS
        + "ref=8 nt=4 q=20 r=63" + LS
        + "ref=9 nt=4 q=20 r=63" + LS
        + "ref=10 nt=4 q=20 r=63" + LS
        + "ref=11 nt=4 q=20 r=63" + LS
        + "ref=12 nt=4 q=20 r=63" + LS
        + "ref=13 nt=4 q=20 r=63" + LS
        + "ref=14 nt=4 q=20 r=63" + LS
        + "ref=15 nt=4 q=20 r=63" + LS
        + "ref=16 nt=4 q=20 r=63" + LS
        + "ref=17 nt=4 q=20 r=63" + LS
        + "ref=18 nt=4 q=20 r=63" + LS
        + "ref=19 nt=4 q=20 r=63" + LS
        + "ref=20 nt=4 q=20 r=63" + LS
        + "ref=21 nt=4 q=20 r=63" + LS
        + "ref=22 nt=4 q=20 r=63" + LS
        + "ref=23 nt=4 q=20 r=63" + LS
        + "ref=24 nt=4 q=20 r=63" + LS
        + "ref=30 nt=1 q=20 r=63" + LS
        + "ref=31 nt=1 q=20 r=63" + LS
        + "ref=32 nt=1 q=20 r=63" + LS
        + "ref=33 nt=1 q=20 r=63" + LS
        + "ref=34 nt=1 q=20 r=63" + LS
        + "ref=35 nt=1 q=20 r=63" + LS
        + "ref=36 nt=1 q=20 r=63" + LS
        + "ref=37 nt=1 q=20 r=63" + LS
        + "ref=38 nt=1 q=20 r=63" + LS
        + "ref=39 nt=1 q=20 r=63" + LS
        , matcher.toString());
  }

  public void testSuperCigarNoTrim() throws Exception {
    //36405 131 paolo-bac 3319  55  10=5N17=1D6=  = 3014  -305  ANCTCCAGGGGTGATCTTCCCACCTCACCTCCC ,!0604)552.665111057,1-.043./5044 AS:i:3  MQ:i:255  XU:Z:1=1R8=5N17=1D2=1X2B5=  XQ:Z:.* XR:Z:A  XA:i:4  RG:Z:CG IH:i:1  NH:i:1

    final byte[] tmpl = DnaUtils.encodeString("GGANCTCCAGGGTTTTTGTGATCTTCCCACCTCAGCCTCCCTTT");
    final SAMFileHeader sft3 = new SAMFileHeader();
    final SAMRecord rec = new SAMRecord(sft3);
    rec.setReadName("blah");
    rec.setReadString("ANCTCCAGGGGTGATCTTCCCACCTCACCTCCC");       //used just for the length
    rec.setBaseQualityString(",!0604)552.665111057,1-.043./5044");
    rec.setReferenceName("sdjr");
    rec.setFlags(131);
    rec.setAlignmentStart(3);
    rec.setInferredInsertSize(0); // Default, overwritten for valid paired-end
    rec.setMateAlignmentStart(0); // Default, overwritten for valid paired-end
    rec.setMateReferenceName("*");
    rec.setCigarString("10=5N17=1D6=");
    rec.setAttribute(SamUtils.CG_SUPER_CIGAR, "1=1R8=5N17=1D2=1X2B5=");
    rec.setAttribute(SamUtils.CG_SUPER_CIGAR_OVERLAP_QUALITY, ".*");
    rec.setAttribute(SamUtils.CG_READ_DELTA, "A");

    //orig read = ANCTCCAGGGGTGATCTTCCCACCTCACCACTCCC
    //orig qual = ,!0604)552.665111057,1-.043..*/5044

    final VariantAlignmentRecord var = new VariantAlignmentRecord(rec);

    final VariantParams vp = new VariantParamsBuilder().ignoreQualityScores(false).create();
    try (MockModelMatcher matcher = new MockModelMatcher()) {
      final CigarParserModel cpm = new CigarParserModel(matcher, matcher, 0, tmpl.length, vp);
      cpm.toMatcher(var, MachineType.COMPLETE_GENOMICS, 20, tmpl);

      assertEquals(""
        + "ref=2 nt=1 q=11 r=0" + LS  //A,
        + "ref=3 nt=0 q=0 r=0" + LS //N!
        + "ref=4 nt=2 q=15 r=0" + LS  //C0
        + "ref=5 nt=4 q=21 r=0" + LS  //T6
        + "ref=6 nt=2 q=15 r=0" + LS  //C0
        + "ref=7 nt=2 q=19 r=0" + LS  //C4
        + "ref=8 nt=1 q=8 r=0" + LS
        + "ref=9 nt=3 q=20 r=0" + LS
        + "ref=10 nt=3 q=20 r=0" + LS
        + "ref=11 nt=3 q=17 r=0" + LS
        + "ref=17 nt=3 q=13 r=0" + LS
        + "ref=18 nt=4 q=21 r=0" + LS
        + "ref=19 nt=3 q=21 r=0" + LS
        + "ref=20 nt=1 q=20 r=0" + LS
        + "ref=21 nt=4 q=16 r=0" + LS
        + "ref=22 nt=2 q=16 r=0" + LS
        + "ref=23 nt=4 q=16 r=0" + LS
        + "ref=24 nt=4 q=15 r=0" + LS
        + "ref=25 nt=2 q=20 r=0" + LS
        + "ref=26 nt=2 q=22 r=0" + LS
        + "ref=27 nt=2 q=11 r=0" + LS
        + "ref=28 nt=1 q=16 r=0" + LS
        + "ref=29 nt=2 q=12 r=0" + LS
        + "ref=30 nt=2 q=13 r=0" + LS
        + "ref=31 nt=4 q=15 r=0" + LS //T
        + "ref=32 nt=2 q=19 r=0" + LS //C
        + "ref=33 nt=1 q=18 r=0" + LS //A
        + "ref=35 nt=2 q=13 r=0" + LS //C.
        + "ref=36 nt=2 q=14 r=0" + LS //C/
        + "ref=37 nt=4 q=20 r=0" + LS //T5
        + (MachineType.CG_TRIM ?  ""
        : "ref=38 nt=2 q=15 r=0" + LS //C0
        + "ref=39 nt=2 q=19 r=0" + LS //C4
        + "ref=40 nt=2 q=19 r=0" + LS), matcher.toString()); //C4

      assertEquals("ref=34 n" + LS, matcher.insertString());
    }
  }


  public void testSuperCigar2() throws Exception {
    //12560913  131 chr21 142 0 10=6N18=1X4=  chr21 13640284  -419  CTCCCAAGTTATTCTCCTGCCTCAGCCTCTCAA 9989::::::<<<<<;=<<;9:;2.*'29/890 AS:i:1  MQ:i:255  XU:Z:10=6N20=2B1X4= XQ:Z:&) XR:Z:C  XA:i:1  RG:Z:GS000005351  IH:i:3  NH:i:3

    final byte[] tmpl = DnaUtils.encodeString("CTCCCAAGTTCAAGTGATTCTCCTGCCTCAGCCTTTCAAGTAGCT");
    final SAMFileHeader sft3 = new SAMFileHeader();
    final SAMRecord rec = new SAMRecord(sft3);
    rec.setReadName("blah");
    rec.setReadString("CTCCCAAGTTATTCTCCTGCCTCAGCCTCTCAA");       //used just for the length
    rec.setBaseQualityString("9989::::::<<<<<;=<<;9:;2.*'29/890");
    rec.setReferenceName("sdjr");
    rec.setFlags(131);
    rec.setAlignmentStart(1);
    rec.setInferredInsertSize(0); // Default, overwritten for valid paired-end
    rec.setMateAlignmentStart(0); // Default, overwritten for valid paired-end
    rec.setMateReferenceName("*");
    rec.setCigarString("10=6N18=1X4=");
    rec.setAttribute(SamUtils.CG_SUPER_CIGAR, "10=6N20=2B1X4=");
    rec.setAttribute(SamUtils.CG_SUPER_CIGAR_OVERLAP_QUALITY, "&)");
    rec.setAttribute(SamUtils.CG_READ_DELTA, "C");

    //orig read = CTCCCAAGTT      ATTCTCCTGCCTCAGCCTTT CTCAA
    //orig qual = 9989::::::      <<<<<;=<<;9:;2.*'2&) 9/890

    final VariantAlignmentRecord var = new VariantAlignmentRecord(rec);

    final VariantParams vp = new VariantParamsBuilder().ignoreQualityScores(false).create();
    try (MockModelMatcher matcher = new MockModelMatcher()) {
      final CigarParserModel cpm = new CigarParserModel(matcher, matcher, 0, tmpl.length, vp);
      cpm.toMatcher(var, MachineType.COMPLETE_GENOMICS, 20, tmpl);

      assertEquals(""
        + "ref=0 nt=2 q=24 r=0" + LS  //C9
        + "ref=1 nt=4 q=24 r=0" + LS  //T9
        + "ref=2 nt=2 q=23 r=0" + LS  //C8
        + "ref=3 nt=2 q=24 r=0" + LS  //C9
        + "ref=4 nt=2 q=25 r=0" + LS
        + "ref=5 nt=1 q=25 r=0" + LS
        + "ref=6 nt=1 q=25 r=0" + LS
        + "ref=7 nt=3 q=25 r=0" + LS
        + "ref=8 nt=4 q=25 r=0" + LS
        + "ref=9 nt=4 q=25 r=0" + LS
        + "ref=16 nt=1 q=27 r=0" + LS
        + "ref=17 nt=4 q=27 r=0" + LS
        + "ref=18 nt=4 q=27 r=0" + LS
        + "ref=19 nt=2 q=27 r=0" + LS
        + "ref=20 nt=4 q=27 r=0" + LS
        + "ref=21 nt=2 q=26 r=0" + LS
        + "ref=22 nt=2 q=28 r=0" + LS
        + "ref=23 nt=4 q=27 r=0" + LS
        + "ref=24 nt=3 q=27 r=0" + LS
        + "ref=25 nt=2 q=26 r=0" + LS
        + "ref=26 nt=2 q=24 r=0" + LS
        + "ref=27 nt=4 q=25 r=0" + LS
        + "ref=28 nt=2 q=26 r=0" + LS
        + "ref=29 nt=1 q=17 r=0" + LS
        + "ref=30 nt=3 q=13 r=0" + LS
        + "ref=31 nt=2 q=9 r=0" + LS
        + "ref=32 nt=2 q=6 r=0" + LS
        + "ref=33 nt=4 q=17 r=0" + LS //T2
        + "ref=34 nt=2 q=24 r=0" + LS //C9
        + "ref=35 nt=4 q=14 r=0" + LS //T/
        + (MachineType.CG_TRIM ?  ""
        : "ref=36 nt=2 q=23 r=0" + LS //C8
        + "ref=37 nt=1 q=24 r=0" + LS //A9
        + "ref=38 nt=1 q=15 r=0" + LS), matcher.toString()); //A0

      assertEquals("", matcher.insertString());
    }
  }

  public void testCgLeftPrune() throws Exception {
//    23667286        67      chr21   9720286 55      20=1X2=6N9=1X   =       9720785 499     TGAGAAACTACTTTGTAATGAGTATCTCACAGT      899798:9777779;+;+',$5(5$8,8:8998       AS:i:2  NM:i:2  MQ:i:255        GS:Z:GAGA       GC:Z:3S2G28S    GQ:Z:90 XA:i:4  RG:Z:GS000005362        IH:i:1  NH:i:1

    final MockModelMatcher matcher = new MockModelMatcher();
    final VariantParams vp = new VariantParamsBuilder().ignoreQualityScores(false).create();
    try {
      final byte[] tmpl = DnaUtils.encodeString("AATGAGAAACTACTTTGTAATGTGTCCCCCCATCTCACAGG");
      final SAMFileHeader sft3 = new SAMFileHeader();
      final SAMRecord rec = new SAMRecord(sft3);
      rec.setReadName("blah");
      rec.setReadString("TGAGAAACTACTTTGTAATGAGTATCTCACAGT");       //used just for the length
      rec.setReferenceName("sdjr");
      rec.setFlags(67);
      rec.setAlignmentStart(3);
      rec.setInferredInsertSize(0); // Default, overwritten for valid paired-end
      rec.setMateAlignmentStart(0); // Default, overwritten for valid paired-end
      rec.setMateReferenceName("*");
      rec.setCigarString("20=1X2=6N9=1X");
      rec.setMappingQuality(63);

      final VariantAlignmentRecord var = new VariantAlignmentRecord(rec);

      new CigarParserModel(matcher, matcher, 0, tmpl.length, vp).toMatcher(var, MachineType.COMPLETE_GENOMICS, 20, tmpl);

    } finally {
      matcher.close();
    }
    assertEquals(""
        + "ref=5 nt=3 q=20 r=63" + LS
        + "ref=6 nt=1 q=20 r=63" + LS
        + "ref=7 nt=1 q=20 r=63" + LS
        + "ref=8 nt=1 q=20 r=63" + LS
        + "ref=9 nt=2 q=20 r=63" + LS
        + "ref=10 nt=4 q=20 r=63" + LS
        + "ref=11 nt=1 q=20 r=63" + LS
        + "ref=12 nt=2 q=20 r=63" + LS
        + "ref=13 nt=4 q=20 r=63" + LS
        + "ref=14 nt=4 q=20 r=63" + LS
        + "ref=15 nt=4 q=20 r=63" + LS
        + "ref=16 nt=3 q=20 r=63" + LS
        + "ref=17 nt=4 q=20 r=63" + LS
        + "ref=18 nt=1 q=20 r=63" + LS
        + "ref=19 nt=1 q=20 r=63" + LS
        + "ref=20 nt=4 q=20 r=63" + LS
        + "ref=21 nt=3 q=20 r=63" + LS
        + "ref=22 nt=1 q=20 r=63" + LS
        + "ref=23 nt=3 q=20 r=63" + LS
        + "ref=24 nt=4 q=20 r=63" + LS
        + "ref=31 nt=1 q=20 r=63" + LS
        + "ref=32 nt=4 q=20 r=63" + LS
        + "ref=33 nt=2 q=20 r=63" + LS
        + "ref=34 nt=4 q=20 r=63" + LS
        + "ref=35 nt=2 q=20 r=63" + LS
        + "ref=36 nt=1 q=20 r=63" + LS
        + "ref=37 nt=2 q=20 r=63" + LS
        + "ref=38 nt=1 q=20 r=63" + LS
        + "ref=39 nt=3 q=20 r=63" + LS
        + "ref=40 nt=4 q=20 r=63" + LS
        , matcher.toString());
  }

  public void testCgLeftPrune2() throws Exception {
//    11619609        179     chr21   9720273 55      8=1X6=1X7=6N10= =       9720764 491     AACAGAAGTTTTCTGGGAAACTAAATGTGTGCA      69:46648:;;;:;<<=<;;<<;9:9:9:99:9       AS:i:2  NM:i:2  MQ:i:255        GS:Z:AGAG       GC:Z:3S2G28S    GQ:Z:-: XA:i:3  RG:Z:GS000005362        IH:i:1  NH:i:1

    final MockModelMatcher matcher = new MockModelMatcher();
    final VariantParams vp = new VariantParamsBuilder().ignoreQualityScores(false).create();
    try {
      final byte[] tmpl = DnaUtils.encodeString("TTAACAGAAGATTTCTGTGAAACTANNNNNNAATGTGTGCA");
      final SAMFileHeader sft3 = new SAMFileHeader();
      final SAMRecord rec = new SAMRecord(sft3);
      rec.setReadName("blah");
      rec.setReadString("AACAGAAGTTTTCTGGGAAACTAAATGTGTGCA");       //used just for the length
      rec.setReferenceName("sdjr");
      rec.setFlags(179);
      rec.setAlignmentStart(3);
      rec.setInferredInsertSize(0); // Default, overwritten for valid paired-end
      rec.setMateAlignmentStart(0); // Default, overwritten for valid paired-end
      rec.setMateReferenceName("*");
      rec.setCigarString("8=1X6=1X7=6N10=");
      rec.setMappingQuality(63);

      final VariantAlignmentRecord var = new VariantAlignmentRecord(rec);

      new CigarParserModel(matcher, matcher, 0, tmpl.length, vp).toMatcher(var, MachineType.COMPLETE_GENOMICS, 20, tmpl);

    } finally {
      matcher.close();
    }
    assertEquals(""
        + "ref=5 nt=1 q=20 r=63" + LS
        + "ref=6 nt=3 q=20 r=63" + LS
        + "ref=7 nt=1 q=20 r=63" + LS
        + "ref=8 nt=1 q=20 r=63" + LS
        + "ref=9 nt=3 q=20 r=63" + LS
        + "ref=10 nt=4 q=20 r=63" + LS
        + "ref=11 nt=4 q=20 r=63" + LS
        + "ref=12 nt=4 q=20 r=63" + LS
        + "ref=13 nt=4 q=20 r=63" + LS
        + "ref=14 nt=2 q=20 r=63" + LS
        + "ref=15 nt=4 q=20 r=63" + LS
        + "ref=16 nt=3 q=20 r=63" + LS
        + "ref=17 nt=3 q=20 r=63" + LS
        + "ref=18 nt=3 q=20 r=63" + LS
        + "ref=19 nt=1 q=20 r=63" + LS
        + "ref=20 nt=1 q=20 r=63" + LS
        + "ref=21 nt=1 q=20 r=63" + LS
        + "ref=22 nt=2 q=20 r=63" + LS
        + "ref=23 nt=4 q=20 r=63" + LS
        + "ref=24 nt=1 q=20 r=63" + LS
        + "ref=31 nt=1 q=20 r=63" + LS
        + "ref=32 nt=1 q=20 r=63" + LS
        + "ref=33 nt=4 q=20 r=63" + LS
        + "ref=34 nt=3 q=20 r=63" + LS
        + "ref=35 nt=4 q=20 r=63" + LS
        + "ref=36 nt=3 q=20 r=63" + LS
        + "ref=37 nt=4 q=20 r=63" + LS
        + "ref=38 nt=3 q=20 r=63" + LS
        + "ref=39 nt=2 q=20 r=63" + LS
        + "ref=40 nt=1 q=20 r=63" + LS
        , matcher.toString());
  }

  public void testCgLeftPruneWithSmallGap() throws Exception {
    final MockModelMatcher matcher = new MockModelMatcher();
    final VariantParams vp = new VariantParamsBuilder().ignoreQualityScores(false).create();
    try {
      final byte[] tmpl = DnaUtils.encodeString("AACCCGTTTTTTTTTAAATTTTTTTTTTGGGGGAAAAAAAAAA");
      final SAMFileHeader sft3 = new SAMFileHeader();
      final SAMRecord rec = new SAMRecord(sft3);
      rec.setReadName("blah");
      rec.setReadString("CCCGTTTTTTTTTTTTTTTTTTTAAAAAAAAAA");       //used just for the length
      rec.setReferenceName("sdjr");
      rec.setFlags(97);
      rec.setAlignmentStart(3);
      rec.setInferredInsertSize(0); // Default, overwritten for valid paired-end
      rec.setMateAlignmentStart(0); // Default, overwritten for valid paired-end
      rec.setMateReferenceName("*");
      rec.setCigarString("13=3N10=5N10=");
      rec.setMappingQuality(63);

      final VariantAlignmentRecord var = new VariantAlignmentRecord(rec);

      new CigarParserModel(matcher, matcher, 0, tmpl.length, vp).toMatcher(var, MachineType.COMPLETE_GENOMICS, 20, tmpl);
    } finally {
      matcher.close();
    }
    // System.out.println(matcher.toString());
    assertEquals(""
        + "ref=5 nt=3 q=20 r=63" + LS
        + "ref=6 nt=4 q=20 r=63" + LS
        + "ref=7 nt=4 q=20 r=63" + LS
        + "ref=8 nt=4 q=20 r=63" + LS
        + "ref=9 nt=4 q=20 r=63" + LS
        + "ref=10 nt=4 q=20 r=63" + LS
        + "ref=11 nt=4 q=20 r=63" + LS
        + "ref=12 nt=4 q=20 r=63" + LS
        + "ref=13 nt=4 q=20 r=63" + LS
        + "ref=14 nt=4 q=20 r=63" + LS
        + "ref=18 nt=4 q=20 r=63" + LS
        + "ref=19 nt=4 q=20 r=63" + LS
        + "ref=20 nt=4 q=20 r=63" + LS
        + "ref=21 nt=4 q=20 r=63" + LS
        + "ref=22 nt=4 q=20 r=63" + LS
        + "ref=23 nt=4 q=20 r=63" + LS
        + "ref=24 nt=4 q=20 r=63" + LS
        + "ref=25 nt=4 q=20 r=63" + LS
        + "ref=26 nt=4 q=20 r=63" + LS
        + "ref=27 nt=4 q=20 r=63" + LS
        + "ref=33 nt=1 q=20 r=63" + LS
        + "ref=34 nt=1 q=20 r=63" + LS
        + "ref=35 nt=1 q=20 r=63" + LS
        + "ref=36 nt=1 q=20 r=63" + LS
        + "ref=37 nt=1 q=20 r=63" + LS
        + "ref=38 nt=1 q=20 r=63" + LS
        + "ref=39 nt=1 q=20 r=63" + LS
        + "ref=40 nt=1 q=20 r=63" + LS
        + "ref=41 nt=1 q=20 r=63" + LS
        + "ref=42 nt=1 q=20 r=63" + LS
        ,  matcher.toString());
  }

  public void testCgRightPrune() throws Exception {
//    43906180        115     chr21   9720271 55      10=5N2=1X20=    =       9719808 -463    TAAACAGAAGTGGGAAACTACTTTGTAATGTGT   9:98:884,601(#525<;8777771),20988       AS:i:1  NM:i:1  MQ:i:255        GS:Z:TGTG       GC:Z:28S2G3S    GQ:Z:9: XA:i:1  RG:Z:GS000005362        IH:i:1  NH:i:1

    final MockModelMatcher matcher = new MockModelMatcher();
    final VariantParams vp = new VariantParamsBuilder().ignoreQualityScores(false).create();
    try {
      final byte[] tmpl = DnaUtils.encodeString("AATAAACAGAAGAAAAATGTGAAACTACTTTGTAATGTGT");
      final SAMFileHeader sft3 = new SAMFileHeader();
      final SAMRecord rec = new SAMRecord(sft3);
      rec.setReadName("blah");
      rec.setReadString("TAAACAGAAGTGGGAAACTACTTTGTAATGTGT");       //used just for the length
      rec.setReferenceName("sdjr");
      rec.setFlags(115);
      rec.setAlignmentStart(3);
      rec.setInferredInsertSize(0); // Default, overwritten for valid paired-end
      rec.setMateAlignmentStart(0); // Default, overwritten for valid paired-end
      rec.setMateReferenceName("*");
      rec.setCigarString("10=5N2=1X20=");
      rec.setMappingQuality(63);

      final VariantAlignmentRecord var = new VariantAlignmentRecord(rec);

      new CigarParserModel(matcher, matcher, 0, tmpl.length, vp).toMatcher(var, MachineType.COMPLETE_GENOMICS, 20, tmpl);

    } finally {
      matcher.close();
    }
    assertEquals(""
        + "ref=2 nt=4 q=20 r=63" + LS
        + "ref=3 nt=1 q=20 r=63" + LS
        + "ref=4 nt=1 q=20 r=63" + LS
        + "ref=5 nt=1 q=20 r=63" + LS
        + "ref=6 nt=2 q=20 r=63" + LS
        + "ref=7 nt=1 q=20 r=63" + LS
        + "ref=8 nt=3 q=20 r=63" + LS
        + "ref=9 nt=1 q=20 r=63" + LS
        + "ref=10 nt=1 q=20 r=63" + LS
        + "ref=11 nt=3 q=20 r=63" + LS
        + "ref=17 nt=4 q=20 r=63" + LS
        + "ref=18 nt=3 q=20 r=63" + LS
        + "ref=19 nt=3 q=20 r=63" + LS
        + "ref=20 nt=3 q=20 r=63" + LS
        + "ref=21 nt=1 q=20 r=63" + LS
        + "ref=22 nt=1 q=20 r=63" + LS
        + "ref=23 nt=1 q=20 r=63" + LS
        + "ref=24 nt=2 q=20 r=63" + LS
        + "ref=25 nt=4 q=20 r=63" + LS
        + "ref=26 nt=1 q=20 r=63" + LS
        + "ref=27 nt=2 q=20 r=63" + LS
        + "ref=28 nt=4 q=20 r=63" + LS
        + "ref=29 nt=4 q=20 r=63" + LS
        + "ref=30 nt=4 q=20 r=63" + LS
        + "ref=31 nt=3 q=20 r=63" + LS
        + "ref=32 nt=4 q=20 r=63" + LS
        + "ref=33 nt=1 q=20 r=63" + LS
        + "ref=34 nt=1 q=20 r=63" + LS
        + "ref=35 nt=4 q=20 r=63" + LS
        + "ref=36 nt=3 q=20 r=63" + LS
        , matcher.toString());
  }

  public void testCgRightPrune2() throws Exception {
//    49451311        131     chr21   9720271 34      10=5N2=1X7=1X12=        =       9719862 -409    TAAACAGAAGTGGGAAACTAATTTGTAATGTGT       *799:,81899<<<;<;;<;::;;;586478:8       AS:i:2  NM:i:2  MQ:i:255        GS:Z:TGTG       GC:Z:28S2G3S    GQ:Z:2: XA:i:3  RG:Z:GS000005362        IH:i:1  NH:i:1

    final MockModelMatcher matcher = new MockModelMatcher();
    final VariantParams vp = new VariantParamsBuilder().ignoreQualityScores(false).create();
    try {
      final byte[] tmpl = DnaUtils.encodeString("AATAAACAGAAGAAAAATGTGAAACTAGTTTGTAATGTGT");
      final SAMFileHeader sft3 = new SAMFileHeader();
      final SAMRecord rec = new SAMRecord(sft3);
      rec.setReadName("blah");
      rec.setReadString("TAAACAGAAGTGGGAAACTAATTTGTAATGTGT");       //used just for the length
      rec.setReferenceName("sdjr");
      rec.setFlags(131);
      rec.setAlignmentStart(3);
      rec.setInferredInsertSize(0); // Default, overwritten for valid paired-end
      rec.setMateAlignmentStart(0); // Default, overwritten for valid paired-end
      rec.setMateReferenceName("*");
      rec.setCigarString("10=5N2=1X7=1X12=");
      rec.setMappingQuality(63);

      final VariantAlignmentRecord var = new VariantAlignmentRecord(rec);

      new CigarParserModel(matcher, matcher, 0, tmpl.length, vp).toMatcher(var, MachineType.COMPLETE_GENOMICS, 20, tmpl);

    } finally {
      matcher.close();
    }
    assertEquals(""
        + "ref=2 nt=4 q=20 r=63" + LS
        + "ref=3 nt=1 q=20 r=63" + LS
        + "ref=4 nt=1 q=20 r=63" + LS
        + "ref=5 nt=1 q=20 r=63" + LS
        + "ref=6 nt=2 q=20 r=63" + LS
        + "ref=7 nt=1 q=20 r=63" + LS
        + "ref=8 nt=3 q=20 r=63" + LS
        + "ref=9 nt=1 q=20 r=63" + LS
        + "ref=10 nt=1 q=20 r=63" + LS
        + "ref=11 nt=3 q=20 r=63" + LS
        + "ref=17 nt=4 q=20 r=63" + LS
        + "ref=18 nt=3 q=20 r=63" + LS
        + "ref=19 nt=3 q=20 r=63" + LS
        + "ref=20 nt=3 q=20 r=63" + LS
        + "ref=21 nt=1 q=20 r=63" + LS
        + "ref=22 nt=1 q=20 r=63" + LS
        + "ref=23 nt=1 q=20 r=63" + LS
        + "ref=24 nt=2 q=20 r=63" + LS
        + "ref=25 nt=4 q=20 r=63" + LS
        + "ref=26 nt=1 q=20 r=63" + LS
        + "ref=27 nt=1 q=20 r=63" + LS
        + "ref=28 nt=4 q=20 r=63" + LS
        + "ref=29 nt=4 q=20 r=63" + LS
        + "ref=30 nt=4 q=20 r=63" + LS
        + "ref=31 nt=3 q=20 r=63" + LS
        + "ref=32 nt=4 q=20 r=63" + LS
        + "ref=33 nt=1 q=20 r=63" + LS
        + "ref=34 nt=1 q=20 r=63" + LS
        + "ref=35 nt=4 q=20 r=63" + LS
        + "ref=36 nt=3 q=20 r=63" + LS
        , matcher.toString());
  }

  // test that right end is pruned when there is a small gap as well as a large one
  public void testCgRightPruneWithSmallGap() throws Exception {
    final MockModelMatcher matcher = new MockModelMatcher();
    final VariantParams vp = new VariantParamsBuilder().ignoreQualityScores(false).create();
    try {
      final byte[] tmpl = DnaUtils.encodeString("AATTTTTTTTTTNNNNNTTTTTTTTTTNTTTAAAAAAACCC");
      final SAMFileHeader sft3 = new SAMFileHeader();
      final SAMRecord rec = new SAMRecord(sft3);
      rec.setReadName("blah");
      rec.setReadString("TTTTTTTTTTTTTTTTTTTTTTTAAAAAAACCC");       //used just for the length
      rec.setReferenceName("sdjr");
      rec.setFlags(131);
      rec.setAlignmentStart(3);
      rec.setInferredInsertSize(0); // Default, overwritten for valid paired-end
      rec.setMateAlignmentStart(0); // Default, overwritten for valid paired-end
      rec.setMateReferenceName("*");
      rec.setCigarString("10=5N10=1N13=");
      rec.setMappingQuality(63);

      final VariantAlignmentRecord var = new VariantAlignmentRecord(rec);

      new CigarParserModel(matcher, matcher, 0, tmpl.length, vp).toMatcher(var, MachineType.COMPLETE_GENOMICS, 20, tmpl);


    } finally {
      matcher.close();
    }
    assertEquals(""
        + "ref=2 nt=4 q=20 r=63" + LS
        + "ref=3 nt=4 q=20 r=63" + LS
        + "ref=4 nt=4 q=20 r=63" + LS
        + "ref=5 nt=4 q=20 r=63" + LS
        + "ref=6 nt=4 q=20 r=63" + LS
        + "ref=7 nt=4 q=20 r=63" + LS
        + "ref=8 nt=4 q=20 r=63" + LS
        + "ref=9 nt=4 q=20 r=63" + LS
        + "ref=10 nt=4 q=20 r=63" + LS
        + "ref=11 nt=4 q=20 r=63" + LS
        + "ref=17 nt=4 q=20 r=63" + LS
        + "ref=18 nt=4 q=20 r=63" + LS
        + "ref=19 nt=4 q=20 r=63" + LS
        + "ref=20 nt=4 q=20 r=63" + LS
        + "ref=21 nt=4 q=20 r=63" + LS
        + "ref=22 nt=4 q=20 r=63" + LS
        + "ref=23 nt=4 q=20 r=63" + LS
        + "ref=24 nt=4 q=20 r=63" + LS
        + "ref=25 nt=4 q=20 r=63" + LS
        + "ref=26 nt=4 q=20 r=63" + LS
        + "ref=28 nt=4 q=20 r=63" + LS
        + "ref=29 nt=4 q=20 r=63" + LS
        + "ref=30 nt=4 q=20 r=63" + LS
        + "ref=31 nt=1 q=20 r=63" + LS
        + "ref=32 nt=1 q=20 r=63" + LS
        + "ref=33 nt=1 q=20 r=63" + LS
        + "ref=34 nt=1 q=20 r=63" + LS
        + "ref=35 nt=1 q=20 r=63" + LS
        + "ref=36 nt=1 q=20 r=63" + LS
        + "ref=37 nt=1 q=20 r=63" + LS
        , matcher.toString());
  }

  public void checkCgPrune(String read, String cigar, String template, boolean forward, int count) throws Exception {
    final MockModelMatcher matcher = new MockModelMatcher();
    final VariantParams vp = new VariantParamsBuilder().ignoreQualityScores(false).create();
    try {
      final byte[] tmpl = DnaUtils.encodeString(template);
      final SAMFileHeader sft3 = new SAMFileHeader();
      final SAMRecord rec = new SAMRecord(sft3);
      rec.setReadName("blah");
      rec.setReadString(read);       //used just for the length
      rec.setReferenceName("sdjr");
      rec.setFlags(forward ? 131 : 67);
      rec.setAlignmentStart(3);
      rec.setInferredInsertSize(0); // Default, overwritten for valid paired-end
      rec.setMateAlignmentStart(0); // Default, overwritten for valid paired-end
      rec.setMateReferenceName("*");
      rec.setCigarString(cigar);
      rec.setMappingQuality(63);

      final VariantAlignmentRecord var = new VariantAlignmentRecord(rec);
      new CigarParserModel(matcher, matcher, 0, tmpl.length, vp).toMatcher(var, MachineType.COMPLETE_GENOMICS, 20, tmpl);
    } finally {
      matcher.close();
    }
    assertFalse(matcher.toString().contains("nt=2"));
    assertFalse(matcher.insertString().contains("nt=2"));
    assertEquals(count, matcher.count());
  }

  // In following tests C nucleotides should not appear in pruned cases
  public void testCgPrune() throws Exception {
    checkCgPrune("TTTTTTTATTTTTTTTTTTTTTTAAAAAAACCC", "7=1X2=5N23=", "AATTTTTTTTTTTTTTTTTTTTTTTTTTTTAAAAAAACCC", true, 30);
    checkCgPrune("TTTTTTTATTTTTTTTTTTTTTTAAAAAAAACCC", "7=1X2=5N24=", "AATTTTTTTGTTNNNNNTTTTTTTTTTTTTAAAAAAAACCC", true, 30);
    checkCgPrune("TTTTTTTATTTTTTTTTTTTTTTAAAAAAAAACCC", "7=1X2=5N25=", "AATTTTTTTGTTNNNNNTTTTTTTTTTTTTAAAAAAAAACCC", true, 30);
    checkCgPrune("TTTTTTTGTTTTTTTTTTTTTTTAAAAAAAACCC", "7=1D2=5N24=", "AATTTTTTTGGTTNNNNNTTTTTTTTTTTTTAAAAAAAACCC", true, 30);
    checkCgPrune("TTTTTTTTTTTTTTTTTTTTTTTAAAAAAAACCC", "7=1I2=5N24=", "AATTTTTTTTTNNNNNTTTTTTTTTTTTTAAAAAAAACCC", true, 29);
    checkCgPrune("CCCTTTTTTTTTTTTTTTTTTTTAAAAAAAAAA", "23=6N10=", "AACCCTTTTTTTTTTTTTTTTTTTTNNNNNNAAAAAAAAAA", false, 30);
    checkCgPrune("CCCTTTTTTTTTTTTTTTTTTTTAAAAAAAAAA", "23=7N10=", "AACCCTTTTTTTTTTTTTTTTTTTTNNNNNNNAAAAAAAAAA", false, 30);
    checkCgPrune("CCCGTTTTTTTTTTTTTTTTTTTTAAAAAAAAAA", "24=7N10=", "AACCCGTTTTTTTTTTTTTTTTTTTTNNNNNNNAAAAAAAAAA", false, 30);
    checkCgPrune("CCCGGTTTTTTTTTTTTTTTTTTTTAAAAAAAAAA", "25=7N10=", "AACCCGGTTTTTTTTTTTTTTTTTTTTNNNNNNNAAAAAAAAAA", false, 30);
    checkCgPrune("CCCTTTTTTTTTTTTTTTTTTTTAAAAAAAAAA", "23=7N2I8=", "AACCCTTTTTTTTTTTTTTTTTTTTNNNNNNNAAAAAAAA", false, 28);
    checkCgPrune("CCCTTTTTTTTTTTTTTTTTTTTAAAAAAAAAA", "23=6N5=5S", "AACCCTTTTTTTTTTTTTTTTTTTTNNNNNNAAAAA", false, 25);
    checkCgPrune("TTTTTTTGTTTTTTTTTTTTTTTAAAAAAACCC", "7S1X2=5N23=", "AACTTNNNNNTTTTTTTTTTTTTAAAAAAACCC", true, 23);
  }

  public void testCgTrimOlap1() throws Exception {
//    31386   67      paolo-bac       33604   55      24=6N10=        =       33997   393     GGTGGGCACCGGTGCCTCGCCCCAAGGGCCAAGG        0353,7311-8/$152834,1'4,52-.*14505      AS:i:0  MQ:i:255        XU:Z:5=1B20=6N10=       XQ:Z:/  XA:i:1  IH:i:1  NH:i:1

    final MockModelMatcher matcher = new MockModelMatcher();
    final VariantParams vp = new VariantParamsBuilder().ignoreQualityScores(false).create();
    try {
      final byte[] tmpl = DnaUtils.encodeString("AAGGTGGGCACCGGTGCCTCGCCCCANNNNNNAGGGCCAAGG");
      final SAMFileHeader sft3 = new SAMFileHeader();
      final SAMRecord rec = new SAMRecord(sft3);
      rec.setReadName("blah");
      rec.setReadString("GGTGGGCACCGGTGCCTCGCCCCAAGGGCCAAGG");       //used just for the length
      rec.setReferenceName("sdjr");
      rec.setFlags(67);
      rec.setAlignmentStart(3);
      rec.setInferredInsertSize(0); // Default, overwritten for valid paired-end
      rec.setMateAlignmentStart(0); // Default, overwritten for valid paired-end
      rec.setMateReferenceName("*");
      rec.setCigarString("24=6N10=");
      rec.setMappingQuality(63);
      rec.setAttribute(SamUtils.CG_SUPER_CIGAR, "5=1B20=6N10=");
//      rec.setAttribute(SamUtils.CG_SUPER_CIGAR_OVERLAP_QUALITY, "&)");
//      rec.setAttribute(SamUtils.CG_READ_DELTA, "C");
//      5=1B20=6N10=

      final VariantAlignmentRecord var = new VariantAlignmentRecord(rec);

      new CigarParserModel(matcher, matcher, 0, tmpl.length, vp).toMatcher(var, MachineType.COMPLETE_GENOMICS, 20, tmpl);

    } finally {
      matcher.close();
    }
    assertEquals(""
        + "ref=6 nt=3 q=20 r=63" + LS
        + "ref=7 nt=3 q=20 r=63" + LS
        + "ref=8 nt=2 q=20 r=63" + LS
        + "ref=9 nt=1 q=20 r=63" + LS
        + "ref=10 nt=2 q=20 r=63" + LS
        + "ref=11 nt=2 q=20 r=63" + LS
        + "ref=12 nt=3 q=20 r=63" + LS
        + "ref=13 nt=3 q=20 r=63" + LS
        + "ref=14 nt=4 q=20 r=63" + LS
        + "ref=15 nt=3 q=20 r=63" + LS
        + "ref=16 nt=2 q=20 r=63" + LS
        + "ref=17 nt=2 q=20 r=63" + LS
        + "ref=18 nt=4 q=20 r=63" + LS
        + "ref=19 nt=2 q=20 r=63" + LS
        + "ref=20 nt=3 q=20 r=63" + LS
        + "ref=21 nt=2 q=20 r=63" + LS
        + "ref=22 nt=2 q=20 r=63" + LS
        + "ref=23 nt=2 q=20 r=63" + LS
        + "ref=24 nt=2 q=20 r=63" + LS
        + "ref=25 nt=1 q=20 r=63" + LS
        + "ref=32 nt=1 q=20 r=63" + LS
        + "ref=33 nt=3 q=20 r=63" + LS
        + "ref=34 nt=3 q=20 r=63" + LS
        + "ref=35 nt=3 q=20 r=63" + LS
        + "ref=36 nt=2 q=20 r=63" + LS
        + "ref=37 nt=2 q=20 r=63" + LS
        + "ref=38 nt=1 q=20 r=63" + LS
        + "ref=39 nt=1 q=20 r=63" + LS
        + "ref=40 nt=3 q=20 r=63" + LS
        + "ref=41 nt=3 q=20 r=63" + LS
        , matcher.toString());
  }
}
