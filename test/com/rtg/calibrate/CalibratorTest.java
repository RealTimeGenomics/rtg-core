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

package com.rtg.calibrate;

import static com.rtg.util.StringUtils.LS;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;

import com.rtg.mode.DnaUtils;
import com.rtg.mode.SequenceType;
import com.rtg.reader.FastaUtils;
import com.rtg.reader.MockArraySequencesReader;
import com.rtg.reader.MockSequencesReader;
import com.rtg.sam.BadSuperCigarException;
import com.rtg.sam.SamUtils;
import com.rtg.util.StringUtils;
import com.rtg.util.TestUtils;
import com.rtg.util.cli.CommandLine;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.NoTalkbackSlimException;
import com.rtg.util.intervals.ReferenceRegions;
import com.rtg.util.intervals.RegionRestriction;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.MemoryPrintStream;
import com.rtg.util.test.FileHelper;

import htsjdk.samtools.SAMRecord;
import junit.framework.TestCase;

/**
 */
public class CalibratorTest extends TestCase {
  static final String EXPECTED_CALIBRATION = ""
    + "#CL\tfoo bar" + LS
    + "@nh:group1\t0\t5" + LS
    + "@nh:group2\t0\t2" + LS
    + "@covar\treadgroup\treadposition:7\tbasequality\tsequence\tequal\tdiff\tins\tdel" + LS
    + "group1\t0\t35\tsequence1\t3\t0\t0\t0" + LS
    + "group1\t0\t35\tsequence2\t2\t0\t0\t0" + LS
    + "group1\t1\t0\tsequence1\t3\t0\t0\t0" + LS
    + "group1\t1\t0\tsequence2\t2\t0\t0\t0" + LS
    + "group1\t2\t35\tsequence2\t2\t0\t0\t0" + LS
    + "group1\t3\t32\tsequence2\t2\t0\t0\t0" + LS
    + "group1\t4\t35\tsequence2\t2\t0\t0\t0" + LS
    + "group1\t5\t0\tsequence2\t2\t0\t0\t0" + LS
    + "group1\t6\t0\tsequence2\t2\t0\t0\t0" + LS
    + "group2\t0\t35\tsequence1\t2\t0\t0\t0" + LS
    + "group2\t1\t0\tsequence1\t2\t0\t0\t0" + LS
    + "group2\t2\t35\tsequence1\t2\t0\t0\t0" + LS
    + "group2\t3\t32\tsequence1\t2\t0\t0\t0" + LS
    + "group2\t4\t35\tsequence1\t2\t0\t0\t0" + LS
    ;


  private File mDir;
  private long mEqualCount;
  private long mDiffCount;

  @Override
  public void setUp() throws Exception {
    Diagnostic.setLogStream();
    CommandLine.setCommandArgs("foo", "bar");
    mDir = FileUtils.createTempDir("test", "calibrator");
  }

  @Override
  public void tearDown() {
    CommandLine.clearCommandArgs();
    assertTrue(FileHelper.deleteAll(mDir));
    mDir = null;
  }

  public static String stripVersion(String calibrateFile) {
    assert calibrateFile.startsWith("#Version");
    return calibrateFile.substring(calibrateFile.indexOf("#CL"));
  }

  public void testToString() {
    final Calibrator cal0 = new Calibrator(new Covariate[] {}, null);
    assertEquals("", cal0.toString());

    final Calibrator cal2 = new Calibrator(new Covariate[] {new CovariateBaseQuality(), new CovariateReadPos(35)}, null);
    assertEquals("basequality\treadposition:35", cal2.toString());
  }

  public void testGetCovariateIndex() {
    final Calibrator cal = new Calibrator(new Covariate[] {new CovariateBaseQuality(), new CovariateReadPos(35)}, null);
    assertEquals(-1, cal.getCovariateIndex(CovariateEnum.READGROUP));
    assertEquals(0, cal.getCovariateIndex(CovariateEnum.BASEQUALITY));
    assertEquals(1, cal.getCovariateIndex(CovariateEnum.READPOSITION));
  }

  public void testProcessStats() throws IOException {
    final Calibrator cal = new Calibrator(new Covariate[] {new CovariateBaseQuality(), new CovariateReadPos(35)}, null);
    add2Records(cal);
    final Calibrator.QuerySpec query = cal.new QuerySpec();
    mEqualCount = mDiffCount = 0;
    final StatsProcessor proc = new StatsProcessor() {
      @Override
      public void process(int[] covariateValues, CalibrationStats stats) {
        if (stats != null) {
          mEqualCount += stats.getEqual();
          mDiffCount += stats.getDifferent();
        }
      }
    };
    cal.processStats(proc, query);
    assertEquals(57, mEqualCount);
    assertEquals(10, mDiffCount);
    CalibrationStats stats = cal.getSums(CovariateEnum.READGROUP, "foo");
    assertEquals(57, stats.getEqual());

    // now do a more restricted query
    assertTrue(query.setValue(CovariateEnum.READPOSITION, 10));
    mEqualCount = mDiffCount = 0;
    cal.processStats(proc, query);
    assertEquals(2, mEqualCount);
    assertEquals(0, mDiffCount);
    stats = cal.getSums(CovariateEnum.READPOSITION, "10");
    assertEquals(2, stats.getEqual());

    // now do a very restricted query
    assertTrue(query.setValue(CovariateEnum.BASEQUALITY, 16));
    mEqualCount = mDiffCount = 0;
    cal.processStats(proc, query);
    assertEquals(1, mEqualCount);
    assertEquals(0, mDiffCount);

    // check non existent returns false
    assertFalse(query.setValue(CovariateEnum.READGROUP, 0));
  }

//  public void testHasReadgroup() {
//    final Calibrator cal = new Calibrator(new Covariate[] {new CovariateBaseQuality()});
//    final MockSamBamRecord sam = new MockSamBamRecord();
//    sam.setField(SamBamConstants.CIGAR_FIELD, "35=");
//    sam.setField(SamBamConstants.SEQ_FIELD, "acgtacgtggacgtacgtggacgtacgtggttttt");
//    sam.setField(SamBamConstants.QUAL_FIELD, "012345678901234567890123456789ABCDE");
//    sam.setField(SamBamConstants.POS_FIELD, "1");
//    sam.setField(SamBamConstants.MAPQ_FIELD, "1");
//    try {
//      cal.processRead(sam);
//      fail();
//    } catch (final NoTalkbackSlimException ex) {
//      assertEquals("quality calibration requires a read group to be specified", ex.getMessage());
//    }
//  }

  public void testHasReadgroup2() {
    final Calibrator cal = new Calibrator(new Covariate[] {new CovariateBaseQuality()}, null);
    final SAMRecord sam = new SAMRecord(null);
    sam.setCigarString("35=");
    sam.setReadString("acgtacgtggacgtacgtggacgtacgtggttttt");
    sam.setBaseQualityString("012345678901234567890123456789ABCDE");
    sam.setAlignmentStart(1);
    sam.setMappingQuality(1);
    try {
      cal.processRead(sam);
      fail();
    } catch (final NoTalkbackSlimException ex) {
      assertEquals("quality calibration requires a read group to be specified", ex.getMessage());
    }
  }

  public void testWriteToFile() throws IOException {
    final Calibrator cal = new Calibrator(new Covariate[] {new CovariateBaseQuality(), new CovariateReadPos(35)}, null);
    add2Records(cal);
    // now check the output file
    final File tmp = File.createTempFile("test", "calibrator", mDir);
    cal.writeToFile(tmp);
    final String stats = FileUtils.fileToString(tmp);
    //System.out.println(stats);
    check2Records(stats);
  }

  public void testWriteToFileLegacy() throws IOException {
    final Calibrator cal = new Calibrator(new Covariate[] {new CovariateBaseQuality(), new CovariateReadPos(35)}, null);
    add2RecordsLegacy(cal);
    // now check the output file
    final File tmp = File.createTempFile("test", "calibrator", mDir);
    cal.writeToFile(tmp);
    final String stats = FileUtils.fileToString(tmp);
//    System.out.println(stats);
    check2Records(stats);
  }

  // Tests that records with bad super cigars are logged and ignored.
  public void testBadCigar() {
    final ByteArrayOutputStream bos = new ByteArrayOutputStream();
    Diagnostic.setLogStream(new PrintStream(bos));
    final Calibrator cal = new Calibrator(new Covariate[] {new CovariateReadPos(3), new CovariateBaseQuality()}, null);
    final SAMRecord sam = new SAMRecord(null);
    sam.setCigarString("3X");
    sam.setAttribute("XU", "3X");
    sam.setAttribute("XR", "TT");
    sam.setReadString("acg");
    sam.setAlignmentStart(1);
    sam.setBaseQualityString("D!D123456789012345678901234567890DD");
    sam.setMappingQuality(1);
    sam.setFlags(67);
    sam.setAttribute("RG", "rg1");
    cal.processRead(sam);
    Diagnostic.setLogStream(); // closes the bos log stream.
    assertTrue(bos.toString(), bos.toString().contains("Ignored SAM record due to Ill-formed cigar/read delta: 3X/TT"));
  }

  // Tests the case where the last index is bigger than the first index
  public void testMax() throws IOException {
    final Calibrator cal = new Calibrator(new Covariate[] {new CovariateReadPos(3), new CovariateBaseQuality()}, null);
    final SAMRecord sam = new SAMRecord(null);
    sam.setCigarString("3=");
    sam.setReadString("acg");
    sam.setAlignmentStart(1);
    sam.setBaseQualityString("D!D");
    sam.setMappingQuality(1);
    sam.setAttribute("RG", "read group 2");
    cal.processRead(sam);
    // now check the output file
    final File tmp = File.createTempFile("test", "calibrator1", mDir);
    cal.writeToFile(tmp);
    String stats = stripVersion(FileUtils.fileToString(tmp));
    //System.out.println(stats);
    final String expected = ""
      + "#CL foo bar "
      + "@nh:read group 2 0 1 "
      + "@covar readposition:3 basequality equal diff ins del "
      + "0 35 1 0 0 0 "
      + "1 0 1 0 0 0 "
      + "2 35 1 0 0 0 ";
    assertEquals(expected, stats.replaceAll("\\s+", " "));
    final Calibrator cal0 = new Calibrator(new Covariate[] {new CovariateReadPos(3), new CovariateBaseQuality()}, null);
    cal0.accumulate(tmp);
    final File tmp2 = File.createTempFile("test", "calibrator2", mDir);
    cal0.writeToFile(tmp2);
    stats = stripVersion(FileUtils.fileToString(tmp2));
    //System.out.println(stats);
    assertEquals(expected, stats.replaceAll("\\s+", " "));
  }

//  public void testExpanding() {
//    final Calibrator cal = new Calibrator(new Covariate[] {new CovariateReadGroup(), new CovariateReadPos(2), new CovariateBaseQuality(), new CovariateSequence()});
//    MockSamBamRecord sam = new MockSamBamRecord();
//    sam.setField(SamBamConstants.RNAME_FIELD, "sequence1");
//    sam.setField(SamBamConstants.CIGAR_FIELD, "2=");
//    sam.setField(SamBamConstants.SEQ_FIELD,  "ac");
//    sam.setField(SamBamConstants.QUAL_FIELD, "D!");
//    sam.setField(SamBamConstants.POS_FIELD, "1");
//    sam.setField(SamBamConstants.MAPQ_FIELD, "1");
//    sam.addAttribute("RG:Z:group1");
//    cal.processRead(sam);
//    cal.processRead(sam);
//    cal.processRead(sam);
//    sam = new MockSamBamRecord();
//    sam.setField(SamBamConstants.RNAME_FIELD, "sequence1");
//    sam.setField(SamBamConstants.CIGAR_FIELD, "5=");
//    sam.setField(SamBamConstants.SEQ_FIELD,  "actga");
//    sam.setField(SamBamConstants.QUAL_FIELD, "D!DAD");
//    sam.setField(SamBamConstants.POS_FIELD, "1");
//    sam.setField(SamBamConstants.MAPQ_FIELD, "1");
//    sam.addAttribute("RG:Z:group2");
//    cal.processRead(sam);
//    cal.processRead(sam);
//    sam = new MockSamBamRecord();
//    sam.setField(SamBamConstants.RNAME_FIELD, "sequence2");
//    sam.setField(SamBamConstants.CIGAR_FIELD, "7=");
//    sam.setField(SamBamConstants.SEQ_FIELD,  "actgact");
//    sam.setField(SamBamConstants.POS_FIELD, "1");
//    sam.setField(SamBamConstants.QUAL_FIELD, "D!DAD!!");
//    sam.setField(SamBamConstants.MAPQ_FIELD, "1");
//    sam.addAttribute("RG:Z:group1");
//    cal.processRead(sam);
//    cal.processRead(sam);
//    // now check the output file
//    final File tmp = File.createTempFile("test", "calibrator1", mDir);
//    cal.writeToFile(tmp);
//    String stats = stripVersion(FileUtils.fileToString(tmp));
//    final String expected = EXPECTED_CALIBRATION;
//    // System.out.println(stats);
//    assertEquals(expected, stats);
//    final Calibrator cal0 = new Calibrator(new Covariate[] {new CovariateReadGroup(), new CovariateReadPos(1), new CovariateBaseQuality(), new CovariateSequence()});
//    cal0.accumulate(tmp);
//    final File tmp2 = File.createTempFile("test", "calibrator2", mDir);
//    cal0.writeToFile(tmp2);
//    stats = stripVersion(FileUtils.fileToString(tmp2));
//    //System.out.println(stats);
//    assertEquals(expected, stats);
//  }

  public void testExpanding2() throws IOException {
    final Calibrator cal = new Calibrator(new Covariate[] {new CovariateReadGroup(), new CovariateReadPos(2), new CovariateBaseQuality(), new CovariateSequence()}, null);
    SAMRecord sam = new SAMRecord(null);
    sam.setReferenceName("sequence1");
    sam.setCigarString("2=");
    sam.setReadString("ac");
    sam.setBaseQualityString("D!");
    sam.setAlignmentStart(1);
    sam.setMappingQuality(1);
    sam.setAttribute("RG", "group1");
    cal.processRead(sam);
    cal.processRead(sam);
    cal.processRead(sam);
    sam = new SAMRecord(null);
    sam.setReferenceName("sequence1");
    sam.setCigarString("5=");
    sam.setReadString("actga");
    sam.setBaseQualityString("D!DAD");
    sam.setAlignmentStart(1);
    sam.setMappingQuality(1);
    sam.setAttribute("RG", "group2");
    cal.processRead(sam);
    cal.processRead(sam);
    sam = new SAMRecord(null);
    sam.setReferenceName("sequence2");
    sam.setCigarString("7=");
    sam.setReadString("actgact");
    sam.setBaseQualityString("D!DAD!!");
    sam.setAlignmentStart(1);
    sam.setMappingQuality(1);
    sam.setAttribute("RG", "group1");
    cal.processRead(sam);
    cal.processRead(sam);
    // now check the output file
    final File tmp = File.createTempFile("test", "calibrator1", mDir);
    cal.writeToFile(tmp);
    String stats = stripVersion(FileUtils.fileToString(tmp));
    final String expected = EXPECTED_CALIBRATION;
    // System.out.println(stats);
    assertEquals(expected, stats);
    final Calibrator cal0 = new Calibrator(new Covariate[] {new CovariateReadGroup(), new CovariateReadPos(1), new CovariateBaseQuality(), new CovariateSequence()}, null);
    cal0.accumulate(tmp);
    final File tmp2 = File.createTempFile("test", "calibrator2", mDir);
    cal0.writeToFile(tmp2);
    stats = stripVersion(FileUtils.fileToString(tmp2));
    //System.out.println(stats);
    assertEquals(expected, stats);
  }

  public void testNotExpanding() throws IOException {
    final Covariate cov = new CovariateImpl("Dummy", 1) {
      @Override
      public int value(SAMRecord sam, CalibratorCigarParser parser) {
        return 0;
      }
      @Override
      public void resized() {
        fail("Should not have resized");
      }

      @Override
      public CovariateEnum getType() {
        return null;
      }

    };
    final Calibrator cal = new Calibrator(new Covariate[] {cov}, null);
    final SAMRecord sam = new SAMRecord(null);
    sam.setCigarString("2=");
    sam.setReadString("ac");
    sam.setAlignmentStart(1);
    sam.setMappingQuality(1);
    sam.setBaseQualityString("*");
    sam.setAttribute("RG", "group1");
    cal.processRead(sam);
    final File tmp = File.createTempFile("test", "calibrator1", mDir);
    cal.writeToFile(tmp);
    cal.accumulate(tmp);
  }

  public void testReadFromFile() throws IOException {
    final Calibrator cal = new Calibrator(new Covariate[] {new CovariateBaseQuality(), new CovariateReadPos(35)}, null);
    add2Records(cal);
    // now check the output file
    final File tmp = File.createTempFile("test", "calibrator1", mDir);
    cal.writeToFile(tmp);
    String stats = FileUtils.fileToString(tmp);
    //System.out.println(stats);
    check2Records(stats);
    final Calibrator cal0 = new Calibrator(new Covariate[] {new CovariateBaseQuality(), new CovariateReadPos(35)}, null);
    cal0.accumulate(tmp);
    final File tmp2 = File.createTempFile("test", "calibrator2", mDir);
    cal0.writeToFile(tmp2);
    stats = FileUtils.fileToString(tmp2);
    //System.out.println(stats);
    check2Records(stats);
  }

  public void testErrors() throws IOException {
    final Calibrator cal = new Calibrator(new Covariate[] {new CovariateReadGroup(), new CovariateReadPos(3)}, null);
    final File f = File.createTempFile("test", "file", mDir);
    String line = "blah";
    FileUtils.stringToFile(line, f);
    try {
      cal.accumulate(f);
    } catch (final NoTalkbackSlimException e) {
      assertEquals("calibration file \"" + f.getPath() + "\" has bad header line: " + line, e.getMessage());
    }
    line = "@covar\tblah";
    FileUtils.stringToFile(line, f);
    try {
      cal.accumulate(f);
    } catch (final NoTalkbackSlimException e) {
      assertEquals("calibration file \"" + f.getPath() + "\" does not have correct stats columns: " + line, e.getMessage());
    }
    line = "@covar\treadgroup\tequal\tdiff\tins\tdel";
    FileUtils.stringToFile(line, f);
    try {
      cal.accumulate(f);
    } catch (final NoTalkbackSlimException e) {
      assertEquals("calibration file \"" + f.getPath() + "\" has mismatching covariates: " + line, e.getMessage());
    }
    line = "@covar\treadposition:3\treadgroup\tequal\tdiff\tins\tdel";
    FileUtils.stringToFile(line, f);
    try {
      cal.accumulate(f);
    } catch (final NoTalkbackSlimException e) {
      assertEquals("calibration file \"" + f.getPath() + "\" contains unexpected covariate: readposition", e.getMessage());
    }
    line = "@covar\treadgroup\treadposition:3\tequal\tdiff\tins\tdel\nblah";
    FileUtils.stringToFile(line, f);
    try {
      cal.accumulate(f);
    } catch (final NoTalkbackSlimException e) {
      assertEquals("calibration file \"" + f.getPath() + "\" contains invalid line: blah", e.getMessage());
    }
  }

  /**
   * @param stats expected output after doing add2Records.
   */
  public void check2Records(final String stats) {
    TestUtils.containsAll(stats,
        "#CL\tfoo bar",
        "@nh:rg1\t0\t2",
        "@mnp:rg1\t0\t0\t0\t0\t0\t2",
        "@ins:rg1\t0\t0\t1",
        "@del:rg1\t0\t0\t0\t1",
        "@covar\tbasequality\treadposition:35\tequal\tdiff\tins\tdel",
        "15   0  2 0 0 0".replaceAll("  *", "\t"),
        "15  10  1 0 0 0".replaceAll("  *", "\t"),
        "16  10  1 0 0 0".replaceAll("  *", "\t"),
        "15  20  2 0 0 0".replaceAll("  *", "\t"),
        "16   1  2 0 0 0".replaceAll("  *", "\t"),
        "16  11  2 0 0 0".replaceAll("  *", "\t"),
        "16  21  2 0 0 0".replaceAll("  *", "\t"),
        "17   2  2 0 0 0".replaceAll("  *", "\t"),
        "17  12  2 0 0 0".replaceAll("  *", "\t"),
        "17  22  2 0 0 3".replaceAll("  *", "\t"),
        "18   3  2 0 0 0".replaceAll("  *", "\t"),
        "18  13  2 0 0 0".replaceAll("  *", "\t"),
        "18  23  2 0 0 0".replaceAll("  *", "\t"),
        "19   4  2 0 0 0".replaceAll("  *", "\t"),
        "19  14  2 0 0 0".replaceAll("  *", "\t"),
        "19  24  2 0 0 0".replaceAll("  *", "\t"),
        "20   5  1 1 0 0".replaceAll("  *", "\t"),
        "20  15  1 0 1 0".replaceAll("  *", "\t"),
        "20  25  2 0 0 0".replaceAll("  *", "\t"),
        "21   6  1 1 0 0".replaceAll("  *", "\t"),
        "21  16  1 0 1 0".replaceAll("  *", "\t"),
        "21  26  1 1 0 0".replaceAll("  *", "\t"),
        "22   7  1 1 0 0".replaceAll("  *", "\t"),
        "22  17  2 0 0 0".replaceAll("  *", "\t"),
        "22  27  1 1 0 0".replaceAll("  *", "\t"),
        "23   8  1 1 0 0".replaceAll("  *", "\t"),
        "23  18  2 0 0 0".replaceAll("  *", "\t"),
        "23  28  1 1 0 0".replaceAll("  *", "\t"),
        "24   9  1 1 0 0".replaceAll("  *", "\t"),
        "24  19  2 0 0 0".replaceAll("  *", "\t"),
        "24  29  1 1 0 0".replaceAll("  *", "\t"),
        "32  30  1 1 0 0".replaceAll("  *", "\t"),
        "33  31  2 0 0 0".replaceAll("  *", "\t"),
        "34  32  2 0 0 0".replaceAll("  *", "\t"),
        "35  33  2 0 0 0".replaceAll("  *", "\t"),
        "36  34  1 0 0 0".replaceAll("  *", "\t")
        );
  }

  /**
   * Adds a few example records to the given calibrator.
   * @param cal the Calibrator to add records to.
   */
  public void add2Records(final Calibrator cal) {
    final SAMRecord sam = new SAMRecord(null);
    sam.setCigarString("35=");
    sam.setReadString("acgtacgtggacgtacgtggacgtacgtggttttt");
    sam.setBaseQualityString("012345678901234567890123456789ABCDE");
    sam.setAlignmentStart(1);
    sam.setMappingQuality(1);
    sam.setAttribute("RG", "rg1");
    cal.processRead(sam);
    sam.setCigarString("5=5X5=2I5=3D4=5X3=");
    sam.setBaseQualityString("012345678911234567890123456789ABCD");
    sam.setAlignmentStart(1);
    sam.setMappingQuality(1);
    cal.processRead(sam);
  }

  /**
   * Adds a few example records to the given calibrator.
   * @param cal the Calibrator to add records to.
   */
  public void add2RecordsLegacy(final Calibrator cal) {
    final SAMRecord sam = new SAMRecord(null);
    sam.setCigarString("35M");
    sam.setReadString("acgtacgtggacgtacgtggacgtacgtggttttt");
    sam.setBaseQualityString("012345678901234567890123456789ABCDE");
    sam.setAlignmentStart(1);
    sam.setMappingQuality(1);        //     acgtacgtggacgtacgtggacgtacgtggtttttacgtacgtggacgta
                                     //                                                     cgtgg   acgtacgtggtttt
    final byte[] t = DnaUtils.encodeString("acgtacgtggacgtacgtggacgtacgtggtttttacgtaAAAAAacgtatggacNNNgtacAAAAAtttt");
    cal.setTemplate(t, t.length);
    sam.setAttribute("RG", "rg1");
    cal.processRead(sam);
    sam.setReadString("acgtacgtggacgtacgtggacgtacgtggtttt");  //1 shorter to test qualities better
    sam.setCigarString("15M2I5M3D12M");
    sam.setBaseQualityString("012345678911234567890123456789ABCD");
    sam.setAlignmentStart(36);
    sam.setMappingQuality(1);
    cal.processRead(sam);
  }

  private static final String CG_EXPECTED = ""
    + "@cggap:rg1          0     0    10     0     0     0   10".replaceAll("  *", "\t") + LS
    + "@cggap:rg2  0    20     0     0     0    20".replaceAll("  *", "\t") + LS
    + "@cgover:rg1         0     0     0     0    10".replaceAll("  *", "\t") + LS
    + "@cgover:rg2 0     0    20".replaceAll("  *", "\t") + LS
    + "@mnp:rg1            0     0    10".replaceAll("  *", "\t") + LS
    + "@nh:rg1             0    10".replaceAll("  *", "\t") + LS
    + "@nh:rg2     0    20".replaceAll("  *", "\t") + LS
    + "@covar readgroup equal diff ins del".replaceAll("  *", "\t") + LS
    + "rg2 700     0     0     0".replaceAll("  *", "\t") + LS
    + "rg1         330    20     0     0".replaceAll("  *", "\t") + LS
    ;
  public void testCG() throws IOException {
    final Calibrator cal = new Calibrator(new Covariate[] {new CovariateReadGroup()}, null);
    final byte[] t = DnaUtils.encodeString("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA");
    cal.setTemplate(t, t.length);
    final SAMRecord sam1 = new SAMRecord(null);
    sam1.setAttribute(SamUtils.CG_SUPER_CIGAR, "5=2B10=1N10=5N10=");
    sam1.setAlignmentStart(1);
    sam1.setMappingQuality(1);
//    sam1.setBaseQualityString("*");
    sam1.setFlags(67);
    sam1.setAttribute("RG", "rg2");
    final SAMRecord sam2 = new SAMRecord(null);
    sam2.setAttribute(SamUtils.CG_SUPER_CIGAR, "5=4B10=2N10=6N3=2X5=");
    sam2.setAttribute(SamUtils.CG_READ_DELTA, "TT");
    sam2.setAlignmentStart(1);
    sam2.setMappingQuality(1);
//    sam2.setBaseQualityString("*");
    sam2.setFlags(67);
    sam2.setAttribute("RG", "rg1");
    for (int i = 0; i < 10; i++) {
      cal.processRead(sam1);
      cal.processRead(sam1);
      cal.processRead(sam2);
    }
    // now check the output file
    final File tmp = File.createTempFile("test", "calibratorCG", mDir);
    cal.writeToFile(tmp);
    String stats = FileUtils.fileToString(tmp);
    //System.err.println(CG_EXPECTED);
    //System.err.println("---------");
    //System.err.println(stats);1
    assertEquals(CG_EXPECTED, StringUtils.grep(stats, "^[^#]"));
    final Calibrator cal0 = new Calibrator(new Covariate[] {new CovariateReadGroup()}, null);
    cal0.accumulate(tmp);
    cal0.accumulate(tmp);
    final File tmp2 = File.createTempFile("test", "calibratorCG2", mDir);
    cal0.writeToFile(tmp2);
    stats = FileUtils.fileToString(tmp2);
    //System.out.println(stats);
    TestUtils.containsAll(stats,
        "@cggap:rg1  0 0  20  0 0 0  20".replaceAll("  *", "\t"),
        "@cggap:rg2  0 40  0  0 0 40".replaceAll("  *", "\t"),
        "@cgover:rg1 0 0 0  0 20".replaceAll("  *", "\t"),
        "@cgover:rg2 0 0 40".replaceAll("  *", "\t"),
        "@mnp:rg1  0 0 20".replaceAll("  *", "\t"),
        "@nh:rg1 0 20".replaceAll("  *", "\t"),
        "@nh:rg2 0 40".replaceAll("  *", "\t"),
        "rg2 1400  0 0 0".replaceAll("  *", "\t"),
        "rg1 660 40  0 0".replaceAll("  *", "\t")
    );
  }

  public void testUnknownTemplate() throws Exception {

    final Calibrator cal = new Calibrator(new Covariate[] {new CovariateReadGroup()}, null);
    final SAMRecord sam1 = new SAMRecord(null);
    sam1.setAttribute(SamUtils.CG_SUPER_CIGAR, "1T");
    sam1.setAlignmentStart(1);
    sam1.setMappingQuality(1);
    sam1.setAttribute("RG", "rg1");

    cal.mSamRec = sam1;

    cal.mParser.setCigar("1T", "A");
    cal.mParser.parse();

    assertEquals(0, cal.findStats(cal.mParser).getEqual());
  }

  public void testUnknownRead() throws Exception {

    final Calibrator cal = new Calibrator(new Covariate[] {new CovariateReadGroup()}, null);
    final SAMRecord sam1 = new SAMRecord(null);
    sam1.setAttribute(SamUtils.CG_SUPER_CIGAR, "1R");
    sam1.setAlignmentStart(1);
    sam1.setMappingQuality(1);
    sam1.setAttribute("RG", "rg1");

    cal.mSamRec = sam1;

    cal.mParser.setCigar("1R", "");
    cal.mParser.parse();

    assertEquals(0, cal.findStats(cal.mParser).getEqual());
  }

  public void testGetCovariates() throws IOException {
    final File f = FileUtils.stringToFile(EXPECTED_CALIBRATION, new File(mDir, "cal"));
    final Covariate[] c = Calibrator.getCovariateSet(f);
    assertEquals(4, c.length);
    assertEquals("readgroup", c[0].name());
    assertTrue(c[0] instanceof CovariateReadGroup);
    assertEquals("readposition:7", c[1].name());
    assertTrue(c[1] instanceof CovariateReadPos);
    assertEquals("basequality", c[2].name());
    assertTrue(c[2] instanceof CovariateBaseQuality);
    assertEquals("sequence", c[3].name());
    assertTrue(c[3] instanceof CovariateSequence);
  }

  class CalibratorQualityTestParser extends Calibrator {
    private int mPos = 0;
    byte[] mExp = FastaUtils.asciiToRawQuality("4316%%68883-56+141663,2.3----45/.,2");

    public CalibratorQualityTestParser(Covariate[] vars) {
      super(vars, null);
    }

    void setExpectedQuals(byte[] expQual) {
      mExp = expQual;
    }

    void resetPos() {
      mPos = 0;
    }

    @Override
    protected CalibrationStats findStats(CalibratorCigarParser currPos) throws BadSuperCigarException {
      assertEquals("at pos " + mPos, mExp[mPos], currPos.getCurrentQuality());
      mPos++;
      return new CalibrationStats(null);
    }
  }

  public void testSuperCigarParsing() throws Exception {
    final CalibratorQualityTestParser cal = new CalibratorQualityTestParser(new Covariate[] {new CovariateReadGroup()});
    final CalibratorCigarParser parser = new CalibratorCigarParser(cal);

    assertEquals(20, parser.getCurrentQuality());

    final byte[] expQual = FastaUtils.asciiToRawQuality("4316%%68883-56+141663,2.3----45/.,2");
    cal.setExpectedQuals(expQual);
    parser.setQualities(expQual);
    assertEquals(19, parser.getCurrentQuality());

    parser.setTemplate(DnaUtils.encodeString("ACTGACTGACTGACTGACTGACTAGTAGCTAGCTAGTCGATCGCATCGTAGCTAG"));
    parser.setCigar("10=5N25=", null);
    parser.setQualities(FastaUtils.asciiToRawQuality("4316%%68883-56+141663,2.3----45/.,2"));
    parser.parse();

    parser.setCigar("5=2B20=6N10=", null);
    parser.setQualities(FastaUtils.asciiToRawQuality("4316%%68883-56+141663,2.3----45/.,2"));
    cal.resetPos();
    parser.parse();
  }


  public void testQualityHandling() {

    final byte[][] exp = new byte[1][];

    final Calibrator calwtf = new Calibrator(new Covariate[] {new CovariateReadGroup()}, null);

    final CalibratorCigarParser parser = new CalibratorCigarParser(calwtf) {
      @Override
      void setQualities(final byte[] qualities) {
        assertTrue("Expected: " + Arrays.toString(exp[0]) + " but was: " + Arrays.toString(qualities), Arrays.equals(exp[0], qualities));
      }

      @Override
      public void parse() {
      }
    };

    final Calibrator cal = new Calibrator(new Covariate[] {new CovariateReadGroup()}, null) {
      @Override
      CalibratorCigarParser getParser() {
        return parser;
      }
    };

    final byte[] t = DnaUtils.encodeString("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA");
    cal.setTemplate(t, t.length);
//    final MockSamBamRecord sam1 = new MockSamBamRecord();
//    sam1.addAttribute(SamUtils.CG_SUPER_CIGAR + ":Z:5=2B20=6N10=");
//    sam1.setField(SamBamConstants.POS_FIELD, "1");
//    sam1.setField(SamBamConstants.MAPQ_FIELD, "1");
//    sam1.setField(SamBamConstants.QUAL_FIELD, "*");
//    sam1.setField(SamBamConstants.FLAG_FIELD, "67");
//    sam1.addAttribute("RG:Z:rg1");

    final SAMRecord sam = new SAMRecord(null);
    sam.setReadString("acgtacgtggacgtacgtggacgtacgtggttttt");
    sam.setAlignmentStart(1);
    sam.setMappingQuality(1);
    sam.setFlags(67);
    sam.setAttribute("RG", "rg1");
    sam.setAttribute(SamUtils.CG_SUPER_CIGAR, "5=2B20=6N10=");

    exp[0] = null;

    cal.processRead(sam);

    sam.setBaseQualityString("J316%8883-56+141663,2.3----45/.,2");
    sam.setAttribute(SamUtils.CG_OVERLAP_QUALITY, "%6");

    exp[0] = FastaUtils.asciiToRawQuality("J316%%68883-56+141663,2.3----45/.,2");

    cal.processRead(sam);

    final MemoryPrintStream mps = new MemoryPrintStream();
    Diagnostic.setLogStream(mps.printStream());
    try {
      final SAMRecord sam2 = new SAMRecord(null);
      sam2.setAlignmentStart(1);
      sam2.setMappingQuality(1);
      sam2.setFlags(131);
      sam2.setBaseQualityString("2,./54----3.2,366141+65-3888%613J");
      sam2.setAttribute("RG", "rg1");
      sam2.setAttribute(SamUtils.CG_SUPER_CIGAR, "10=6N20=2B5=");
      exp[0] = FastaUtils.asciiToRawQuality("2,./54----3.2,366141+65-3888%6%613J");
      cal.processRead(sam2);

      assertTrue(mps.toString(), mps.toString().contains("Ignored SAM record due to SAM record qualities plus XQ not expected length. Was: 33 expected: 35"));

      mps.reset();
      sam2.setAttribute(SamUtils.CG_OVERLAP_QUALITY, "%6");
      cal.processRead(sam2);

      assertEquals(mps.toString(), 0, mps.toString().length()); //check there's no error output this time

      final SAMRecord sam3 = new SAMRecord(null);
      sam3.setAlignmentStart(1);
      sam3.setMappingQuality(1);
      sam3.setFlags(179);
      sam3.setBaseQualityString("599(:.5,5:0;:2:<96;8:;6::9#::99938");
      sam3.setAttribute("RG", "rg1");
      sam3.setAttribute(SamUtils.CG_SUPER_CIGAR, "5=1B20=5N2=1X5=1X1=");
      sam3.setAttribute(SamUtils.CG_OVERLAP_QUALITY, "4");

      exp[0] = FastaUtils.asciiToRawQuality("599(:4.5,5:0;:2:<96;8:;6::9#::99938");
      cal.processRead(sam3);

      assertEquals(mps.toString(), 0, mps.toString().length()); //check there's no error output at all

      final SAMRecord sam4 = new SAMRecord(null);
      sam4.setAlignmentStart(1);
      sam4.setMappingQuality(1);
      sam4.setFlags(115);
      sam4.setBaseQualityString("77:(9978885/5:0<59;46677746377978");
      sam4.setAttribute(SamUtils.CG_OVERLAP_QUALITY, "61");
      sam4.setAttribute("RG", "rg1");
      sam4.setAttribute(SamUtils.CG_SUPER_CIGAR, "10=5N7=1I12=2B5=");
      exp[0] = FastaUtils.asciiToRawQuality("77:(9978885/5:0<59;4667774636177978");
      cal.processRead(sam4);

      assertEquals(mps.toString(), 0, mps.toString().length()); //check there's no error output at all
    } finally {
      Diagnostic.setLogStream();
    }
  }


  static final String LENGTHS = ""
                                + "#CL\tfoo bar" + LS
                                + "@nh:group1\t0\t5" + LS
                                + "@nh:group2\t0\t2" + LS
                                + "@sequence\t10000\tsequence1" + LS
                                + "@sequence\t1234\tsequence2" + LS
                                + "@covar\treadgroup\treadposition:7\tbasequality\tsequence\tequal\tdiff\tins\tdel" + LS
                                + "group1\t0\t35\tsequence1\t3\t0\t0\t0" + LS
                                + "group1\t0\t35\tsequence2\t2\t0\t0\t0" + LS
                                + "group1\t1\t0\tsequence1\t3\t0\t0\t0" + LS
                                + "group1\t1\t0\tsequence2\t2\t0\t0\t0" + LS
                                + "group1\t2\t35\tsequence2\t2\t0\t0\t0" + LS
                                + "group1\t3\t32\tsequence2\t2\t0\t0\t0" + LS
                                + "group1\t4\t35\tsequence2\t2\t0\t0\t0" + LS
                                + "group1\t5\t0\tsequence2\t2\t0\t0\t0" + LS
                                + "group1\t6\t0\tsequence2\t2\t0\t0\t0" + LS
                                + "group2\t0\t35\tsequence1\t2\t0\t0\t0" + LS
                                + "group2\t1\t0\tsequence1\t2\t0\t0\t0" + LS
                                + "group2\t2\t35\tsequence1\t2\t0\t0\t0" + LS
                                + "group2\t3\t32\tsequence1\t2\t0\t0\t0" + LS
                                + "group2\t4\t35\tsequence1\t2\t0\t0\t0" + LS;
  public void testLengths() throws IOException {
    final Calibrator cal = new Calibrator(new Covariate[] {new CovariateReadGroup(), new CovariateReadPos(2), new CovariateBaseQuality(), new CovariateSequence()}, null);
    cal.accumulate(new ByteArrayInputStream(LENGTHS.getBytes()), "testLengths");
    assertEquals(10000, cal.getSequenceLengths().get("sequence1").intValue());
    assertEquals(1234, cal.getSequenceLengths().get("sequence2").intValue());
    final MemoryPrintStream mps = new MemoryPrintStream();
    cal.writeToStream(mps.printStream());
    assertEquals(LENGTHS, StringUtils.grepMinusV(mps.toString(), "Version"));
    final Map<String, Integer> dummyLengths = new HashMap<>();
    dummyLengths.put("sequence1", 300);
    cal.setSequenceLengths(dummyLengths);
    assertEquals(300, cal.getSequenceLengths().get("sequence1").intValue());
  }
  public void testBed() throws IOException {
    final ReferenceRegions bed = new ReferenceRegions();
    bed.add("sequence1", 11, 20);
    final Calibrator cal = new Calibrator(new Covariate[] {new CovariateReadGroup(), new CovariateSequence()}, bed);
    final Map<String, Integer> dummyLengths = new HashMap<>();
    dummyLengths.put("sequence1", 10);
    cal.setSequenceLengths(dummyLengths);
    final byte[] t = DnaUtils.encodeString("acgtgtgaggacgtacgtggacgtacgtggacgtacgtggtttttacgtagtagattttagaggaggggagaaaaccacacgagacagtgtgtgt");
    cal.setTemplate(t, t.length);

    SAMRecord sam = new SAMRecord(null);
    sam.setReadString("acgtacgtggacgtacgtggacgtacgtggttttt");
    sam.setAlignmentStart(11);
    sam.setMappingQuality(1);
    sam.setReferenceName("sequence1");
    sam.setFlags(67);
    sam.setAttribute("RG", "rg1");
    sam.setAttribute(SamUtils.CG_SUPER_CIGAR, "5=2B20=6N10=");
    sam.setCigarString("36=");
//    sam.setAlignmentEnd(50);
    cal.processRead(sam);
    sam = new SAMRecord(null);
    sam.setReadString("acgtacgtggacgtacgtggacgtacgtggttttt");
    sam.setAlignmentStart(21);
    sam.setMappingQuality(1);
    sam.setReferenceName("sequence1");
    sam.setFlags(67);
    sam.setAttribute("RG", "rg1");
    sam.setCigarString("36=");
//    sam.setAlignmentEnd(50);
    cal.processRead(sam);
    sam = new SAMRecord(null);
    sam.setReadString("acgtgtgaggggg");
    sam.setAlignmentStart(1);
    sam.setMappingQuality(1);
    sam.setReferenceName("sequence1");
    sam.setFlags(67);
    sam.setAttribute("RG", "rg1");
    sam.setCigarString("12=");
//    sam.setAlignmentEnd(50);
    cal.processRead(sam);
    final File tmp2 = File.createTempFile("test", "calibrator2", mDir);
    cal.writeToFile(tmp2);
    final String expected = "#CL\tfoo bar" + LS
                      + "@cgover:rg1\t0\t0\t1" + LS
                      + "@nh:rg1\t0\t2" + LS
                      + "@sequence\t10\tsequence1" + LS
                      + "@covar\treadgroup\tsequence\tequal\tdiff\tins\tdel" + LS
                      + "rg1\tsequence1\t12\t0\t0\t0" + LS;
    final String actual = stripVersion(FileUtils.fileToString(tmp2));
    assertEquals(expected, actual);
  }

  public void checkLengthMap(MockSequencesReader reader) throws IOException {
    final Map<String, Integer> names = Calibrator.getSequenceLengthMap(reader, (RegionRestriction) null);
    assertEquals(10, names.size());
    for (int i = 0; i < names.size(); i++) {
      final String key = "seq" + i;
      assertTrue(names.containsKey(key));
      assertEquals(i, (int) names.get(key));
    }
  }

  public void testSequenceLengthMap() throws IOException {
    final MockSequencesReader reader = new MockArraySequencesReader(SequenceType.DNA, new int[] {0, 1, 2, 3, 4, 5, 6, 7, 8, 9});
    checkLengthMap(reader);
  }


}
