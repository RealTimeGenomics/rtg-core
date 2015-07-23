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
package com.rtg.simulation.cnv;

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;
import java.io.StringReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import com.rtg.reader.PrereadType;
import com.rtg.reader.ReaderTestUtils;
import com.rtg.reader.SdfWriter;
import com.rtg.reader.SequencesReader;
import com.rtg.reader.SequencesReaderFactory;
import com.rtg.simulation.CnvFileChecker;
import com.rtg.util.Constants;
import com.rtg.util.InvalidParamsException;
import com.rtg.util.PortableRandom;
import com.rtg.util.StringUtils;
import com.rtg.util.TestUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.ErrorType;
import com.rtg.util.io.FileUtils;
import com.rtg.util.test.FileHelper;
import com.rtg.util.test.NotRandomRandom;

import junit.framework.TestCase;

public class CnvSimulatorTest extends TestCase {

  public CnvSimulatorTest(String name) {
    super(name);
  }

  private File mDir;
  private CnvPriorParams mPriors;
  private PortableRandom mRandom;
  private static final String TB = "\t";

  SequencesReader mDsrOne = null;
  SequencesReader mDsrTwo = null;
  File mCnvs = null;
  SdfWriter mOutput = null;
  SdfWriter mTwin = null;

  private final class MyCnvSimRandom extends CnvSimulator {

    public MyCnvSimRandom(SequencesReader dsr, SdfWriter output, SdfWriter twin, OutputStream cnvs, PortableRandom random,
        CnvPriorParams priors, double percent, int count) {
      super(dsr, output, twin, cnvs, random, priors, percent, count);

      mBreakPoints = taskaGetBreakPoints(10, 10, -1, -1, -1);
      TestUtils.containsAll("[1, 3, 4, 6, 7, 9]", Arrays.toString(mBreakPoints));
      mRegionSequences = new ArrayList<>();
      TestUtils.containsAll(this.toString(), "breakpoints 6", "1", "3", "4", "6", "7", "9");
    }
  }

  private final class MyCnvSimOneSequence extends CnvSimulator {

    public MyCnvSimOneSequence(SequencesReader dsr, SdfWriter output, SdfWriter twin, OutputStream cnvs, PortableRandom random,
        CnvPriorParams priors, double percent, int count) throws IOException {
      super(dsr, output, twin, cnvs, random, priors, percent, count);

      assertEquals("50.00%", percentString(20, 40));
      assertTrue(Arrays.equals("AAACCC".getBytes(), getDnaArray("AAAACCCC".getBytes(), 1, 6)));
      TestUtils.containsAll(headerLines(), "#Version ",
        ", cnvsim", "#CL");
      assertEquals(0, chooseFromAccumDistribution(0.5, new double[] {0.6, 1.0}));
      assertEquals(2, chooseFromAccumDistribution(0.5, new double[] {0.1, 0.2, 1.0}));

      assertEquals(2, chooseFromAccumDistribution(0.5, priors.powerLengthThresholds()));
      assertEquals(50001, randomPriorsLength(0.5, 4));
      mSetCountCnved = 1;
      mRegionSequences = new ArrayList<>();

      mBreakPoints = taskaGetBreakPoints(200, 50, -1, -1, -1);
      mSequenceLengths = new int[] {50, 50, 50, 50};
      mNumberSequences = 4;
      initialize();

      taskbGenerateCnvRegionsObjects();

      //TODO: enter that test again
//      try {
//        setLengthsAndCount();
//        mSetCountCnved = 20000;
//        selectCnvsAndPrepareCopying();
//        fail();
//      } catch (NoTalkbackSlimException e) {
//        assertTrue(true);
////        TestUtils.containsAll(e.getMessage(), new String[] {"Cannot find enough regions for ",
////          "CNV regions;", " possible."});
//
//      }
      // histograms
      CnvRegion region = new CnvRegion(1, 1, 2);
      addToHistogram(region, mLengthHistogramCnved);
      assertTrue(Arrays.equals(new int[]{0, 1, 0, 0, 0, 0, 0, 0, 0, 0}, mLengthHistogramCnved));
      region = new CnvRegion(1, 1, 20000000);
      addToHistogram(region, mLengthHistogramCnved);
      assertTrue(Arrays.equals(new int[]{0, 1, 0, 0, 0, 0, 0, 0, 1, 0}, mLengthHistogramCnved));
      region = new CnvRegion(1, 1, 200000000);
      addToHistogram(region, mLengthHistogramCnved);
      assertTrue(Arrays.equals(new int[]{0, 1, 0, 0, 0, 0, 0, 0, 1, 1}, mLengthHistogramCnved));
      addToHistogram(region, mLengthHistogramCnved);
      assertTrue(Arrays.equals(new int[]{0, 1, 0, 0, 0, 0, 0, 0, 1, 2}, mLengthHistogramCnved));
      region = new CnvRegion(1, 1, 1000000000);
      addToHistogram(region, mLengthHistogramCnved);

      // build regions where attempts fail
      mBreakPoints = new long[0];
      mSequenceLengths = new int[]{20};
      mRegionDummyList = new ArrayList<>();
      mRegionSequences = new ArrayList<>();
      taskbGenerateCnvRegionsObjects();
      final String sb = regionSequencesToString();
      TestUtils.containsAll(sb, "",
        "seq: 0 cnved: false start: 0 length: 20 del-1: false del-2: false copies: 0",
        "seq: 0 cnved: false start: 20 length: 0 del-1: false del-2: false copies: 0");
      TestUtils.containsAll(this.toString(), "CnvSimulator",
        "No. breakpoints 0",
        "Sequence 0 No. regions 2",
        "priors set");

      // prepare how many are selected
      mTotalSequenceLength = 19;
      mSetCountCnved = Integer.MAX_VALUE;
      taskcSetLengthsAndCount(0);
      assertEquals(2, mSetCnvedSequenceLength);
      // select and prepare for copying
      mNumberSequences = 1;
      taskdSelectCnvs();
      taskePrepareCopying();

      assertEquals(0, mTotalCnvedCount);
      assertEquals(1, mCountsPerSeq.length);
      assertEquals(0, mCountsPerSeq[0]);
      assertTrue(mTotalCnvedLength == 0);
      assertEquals(1001, mTries);

      final File psFile = new File(mDir, "console");
      try (FileOutputStream ps = new FileOutputStream(psFile)) {
        outputHistograms(ps);
        ps.flush();
      }

      final String psfilestr = FileUtils.fileToString(psFile.getPath());
      TestUtils.containsAll(psfilestr, "",
        "Gave up after 1001 attempts",
        "Sequence" + TB + "0" + TB + "length" + TB + "20" + TB
        + "CNV-count" + TB + "0" + TB + "CNV-percent" + TB + "0.00%"
      );

      // build CNV regions
      mBreakPoints = new long[]{3};
      mSequenceLengths = new int[]{20};
      mRegionDummyList = new ArrayList<>();
      mRegionSequences = new ArrayList<>();
      taskbGenerateCnvRegionsObjects();
      final String sb1 = regionSequencesToString();
      TestUtils.containsAll(sb1, "",
        "seq: 0 cnved: false start: 0 length: 3 del-1: false del-2: false copies: 0",
        "seq: 0 cnved: false start: 3 length: 17 del-1: false del-2: false copies: 0",
        "seq: 0 cnved: false start: 20 length: 0 del-1: false del-2: false copies: 0");

      mBreakPoints = new long[]{20, 40, 60, 80, 100, 120, 140, 160, 180};
      mSequenceLengths = new int[]{200};
      mRegionDummyList = new ArrayList<>();
      mRegionSequences = new ArrayList<>();
      taskbGenerateCnvRegionsObjects();
      mSetCountCnved = Integer.MAX_VALUE;
      mSetCnvedSequenceLength = 30;
      taskdSelectCnvs();
      taskePrepareCopying();

      assertEquals(2, mTotalCnvedCount);
      assertEquals(1, mCountsPerSeq.length);
      assertEquals(2, mCountsPerSeq[0]);
      assertEquals(40, mLengthCnvedPerSeq[0]);
      assertTrue(mTotalCnvedLength > 0);

      mBreakPoints = new long[]{20, 40, 60, 80, 100, 120, 140, 160, 180};
      mSequenceLengths = new int[]{200};
      mRegionDummyList = new ArrayList<>();
      mRegionSequences = new ArrayList<>();
      taskbGenerateCnvRegionsObjects();
      mSetCountCnved = 2;
      mSetCnvedSequenceLength = Long.MAX_VALUE;
      taskdSelectCnvs();
      taskePrepareCopying();

      assertEquals(3, mTotalCnvedCount);
      assertEquals(1, mCountsPerSeq.length);
      assertEquals(3, mCountsPerSeq[0]);
      assertEquals(60, mLengthCnvedPerSeq[0]);
      assertTrue(mTotalCnvedLength > 0);

      final String sb2 = regionSequencesToString();
      //System.err.println(sb2.toString());
      TestUtils.containsAll(sb2, "",
        "cnved: true");
    }
  }

  private final class MyCnvSimulatorTaske extends CnvSimulator {

    public MyCnvSimulatorTaske(SequencesReader dsr, SdfWriter output, SdfWriter twin, OutputStream cnvs, PortableRandom random,
        CnvPriorParams priors, double percent, int count) {
      super(dsr, output, twin, cnvs, random, priors, percent, count);

      mRegionSequences = new ArrayList<>();
      final ArrayList<CnvRegion> list = new ArrayList<>();
      final CnvRegion copied = new CnvRegion(0, 0, 10, mPriors);
      copied.mNumCopies = 3;
      list.add(copied);
      list.add(new CnvRegion(0, 10, 20, mPriors));
      list.add(new CnvRegion(0, 20, 0)); // endRegion
      mRegionSequences.add(list);

      taskePrepareCopying();

      final StringBuilder str  = new StringBuilder();
      for (final List<CnvRegion> mRegionSequence : mRegionSequences) {
        for (final CnvRegion currentRegion : mRegionSequence) {
          //System.err.println(currentRegion.toString());
          str.append(currentRegion.toString()).append(StringUtils.LS);
        }
      }
      TestUtils.containsAll(str.toString(), "copies: 3", "num copies added: 2", "num flags added: 2",
        "heterolist: false true");
    }
  }

  public void testMyCnvSimulatorTaske() throws IOException {
    Diagnostic.setLogStream();
    try (OutputStream os = new FileOutputStream(mCnvs)) {
      new MyCnvSimulatorTaske(mDsrTwo, mOutput, mTwin, os, mRandom, mPriors, 10.0, Integer.MAX_VALUE);
    }
  }

  private final class MyCnvSimTwoSequences extends CnvSimulator {

    public MyCnvSimTwoSequences(SequencesReader dsr, SdfWriter output, SdfWriter twin, OutputStream cnvs, PortableRandom random,
        CnvPriorParams priors, double percent, int count) throws IOException {
      super(dsr, output, twin, cnvs, random, priors, percent, count);

      // build CNV regions
      mBreakPoints = new long[]{25};
      mSequenceLengths = new int[]{20, 10};
      mRegionDummyList = new ArrayList<>();
      mRegionSequences = new ArrayList<>();
      initialize();
      taskbGenerateCnvRegionsObjects();
      final String sb1 = regionSequencesToString();
      TestUtils.containsAll(sb1,
        "seq: 0 cnved: false start: 0 length: 20 del-1: false del-2: false copies: 0",
        "seq: 0 cnved: false start: 20 length: 0 del-1: false del-2: false copies: 0",
        "seq: 1 cnved: false start: 0 length: 5 del-1: false del-2: false copies: 0",
        "seq: 1 cnved: false start: 5 length: 5 del-1: false del-2: false copies: 0",
        "seq: 1 cnved: false start: 10 length: 0 del-1: false del-2: false copies: 0");
    }
  }

  private final class MyCnvSimAccumulate extends CnvSimulator {


    public MyCnvSimAccumulate(SequencesReader dsr, SdfWriter output, SdfWriter twin, OutputStream cnvs, PortableRandom random,
        CnvPriorParams priors, double percent, int count) {
      super(dsr, output, twin, cnvs, random, priors, percent, count);
      mBreakPoints = new long[] {3, 20, 40};
      mSequenceLengths = new int[]{100};
      mNumberSequences = 1;
      mRegionDummyList = new ArrayList<>();
      mRegionSequences = new ArrayList<>();
      taskbGenerateCnvRegionsObjects();
      taskgAccumulateNonCnvRegions();
      final String sb = regionSequencesToString();
      TestUtils.containsAll(sb, "seq: 0 cnved: false start: 0 length: 100 del-1: false del-2: false copies: 0",
        "seq: 0 cnved: false start: 100 length: 0 del-1: false del-2: false copies: 0");
    }
  }


  public void testMyCnvSimOneSequence() throws IOException {
    Diagnostic.setLogStream();
    try (OutputStream os = new FileOutputStream(mCnvs)) {
      new MyCnvSimOneSequence(mDsrOne, mOutput, mTwin, os, mRandom, mPriors, 10.0, Integer.MAX_VALUE);
    }
  } //TODO:

  public void testMyCnvSimTwoSequences() throws IOException {
    Diagnostic.setLogStream();
    try (OutputStream os = new FileOutputStream(mCnvs)) {
      new MyCnvSimTwoSequences(mDsrTwo, mOutput, mTwin, os, mRandom, mPriors, 10.0, Integer.MAX_VALUE);
    }
  }


  public void testMyCnvSimAccumulate() throws IOException {
    Diagnostic.setLogStream();
    try (OutputStream os = new FileOutputStream(mCnvs)) {
      new MyCnvSimAccumulate(mDsrOne, mOutput, mTwin, os, mRandom, mPriors, 10.0, Integer.MAX_VALUE);
    }
  }

  public void testMyCnvSimRandom() throws IOException {
    Diagnostic.setLogStream();
    final NotRandomRandom random = new NotRandomRandom();
    try (OutputStream os = new FileOutputStream(mCnvs)) {
      new MyCnvSimRandom(mDsrOne, mOutput, mTwin, os, random, mPriors, 10.0, Integer.MAX_VALUE);
    }
  }

  @Override
  public void setUp() throws IOException, InvalidParamsException {
    mDir = FileHelper.createTempDirectory();
    mRandom = new NotRandomRandom();
    mPriors = CnvPriorParams.builder().cnvpriors("cnv-default").create();

    final String genomeOne = ""
      + ">a\n"
      + "AAAAATTTTTCCCCCGGGGG" + StringUtils.LS
      ;
    final File inOne = ReaderTestUtils.getDNASubDir(genomeOne, mDir);
    final String genomeTwo = ""
      + ">a\n"
      + "AAAAATTTTTCCCCCGGGGG" + StringUtils.LS
      + ">b\n"
      + "TTTTTCCCCC" + StringUtils.LS;
    final File inTwo = ReaderTestUtils.getDNASubDir(genomeTwo, mDir);
    mDsrOne = SequencesReaderFactory.createDefaultSequencesReader(inOne);
    mDsrTwo = SequencesReaderFactory.createDefaultSequencesReader(inTwo);

    mCnvs = new File(mDir, "test.cnv");
    final File outputDirectory = new File(mDir, "out");
    final File twinDirectory = new File(mDir, "twin");
    mOutput = new SdfWriter(outputDirectory, Constants.MAX_FILE_SIZE,
        PrereadType.UNKNOWN, false, true, false, mDsrOne.type());
    mTwin = new SdfWriter(twinDirectory, Constants.MAX_FILE_SIZE,
        PrereadType.UNKNOWN, false, true, false, mDsrOne.type());
  }

  public void testFixedRegion() {
    final ByteArrayOutputStream bos = new ByteArrayOutputStream();
    final PrintStream ps = new PrintStream(bos);
    Diagnostic.setLogStream(ps);
    try {
      new CnvSimulator.FixedRegion("blah", 0, 1, 0);
      fail();
    } catch (InvalidParamsException e) {
      assertEquals(ErrorType.INVALID_PARAMETER, e.getErrorType());
      assertTrue(e.getMessage().contains("CNV invalid start end positions, start: 0 end: 1"));
    }
    try {
      new CnvSimulator.FixedRegion("blah", 2, 1, 0);
      fail();
    } catch (InvalidParamsException e) {
      assertEquals(ErrorType.INVALID_PARAMETER, e.getErrorType());
      assertTrue(e.getMessage().contains("CNV invalid start end positions, start: 2 end: 1"));
    }
    final CnvSimulator.FixedRegion reg1 = new CnvSimulator.FixedRegion("reg1", 1, 4, 1);
    assertEquals("reg1", reg1.mSequenceName);
    assertEquals(0, reg1.mStartPosNullBased);
    assertEquals(3, reg1.mEndPosNullBased);
    assertEquals(1, reg1.mCopyNumber);
    assertEquals(4, reg1.mLength);
    final CnvSimulator.FixedRegion reg2 = new CnvSimulator.FixedRegion("reg2", 2, 6, 2);
    assertEquals("reg2", reg2.mSequenceName);
    assertEquals(1, reg2.mStartPosNullBased);
    assertEquals(5, reg2.mEndPosNullBased);
    assertEquals(2, reg2.mCopyNumber);
    assertEquals(5, reg2.mLength);
    assertEquals(Integer.valueOf(0).compareTo(1), CnvSimulator.getFixedRegionComparator().compare(reg1, reg2));
    Diagnostic.setLogStream();
  }

  @Override
  public void tearDown() throws IOException  {
    if (mTwin != null) {
      mTwin.close();
    }
    if (mOutput != null) {
      mOutput.close();
    }
    if (mDsrOne != null) {
      mDsrOne.close();
    }
    if (mDsrTwo != null) {
      mDsrTwo.close();
    }
   assertTrue(FileHelper.deleteAll(mDir));
    mDir = null;
    mTwin = null;
    mOutput = null;
    mDsrOne = null;
    mDsrTwo = null;
    mRandom = null;
    mPriors = null;
    mCnvs = null;
  }

  public void testGenerateLength() throws IOException {
    Diagnostic.setLogStream();
    //final File temp = FileHelper.createTempDirectory();
    final String genome = ""
      + ">a\n"
      + "ACGTACGATCAGCATCTGACATGCTAACGGTCATC" + StringUtils.LS;
    final File in = ReaderTestUtils.getDNASubDir(genome, mDir);
    try {
      final SequencesReader dsr = SequencesReaderFactory.createDefaultSequencesReader(in);
      try {
        final File cnvs = new File(mDir, "test.cnv");
        final File outputDirectory = new File(mDir, "out");
        final File twinDirectory = new File(mDir, "twin");

        final SdfWriter output = new SdfWriter(outputDirectory, Constants.MAX_FILE_SIZE,
            PrereadType.UNKNOWN, false, true, false, dsr.type());
        try {
          final SdfWriter twin = new SdfWriter(twinDirectory, Constants.MAX_FILE_SIZE,
              PrereadType.UNKNOWN, false, true, false, dsr.type());
          try {
            final CnvSimulator cs = new CnvSimulator(dsr, output, twin, new FileOutputStream(cnvs), mRandom, mPriors, 10, Integer.MAX_VALUE);
            cs.generate();
            final String csstr = cs.toString();
            //System.err.println(csstr);
            TestUtils.containsAll(csstr, "CnvSimulator", "No. breakpoints ",
              "Sequence 0 No. regions 2", "priors set");
            final File psFile = new File(mDir, "console");
            final FileOutputStream ps = new FileOutputStream(psFile);
            try {
              cs.outputHistograms(ps);
            } finally {
              ps.flush();
              ps.close();
              final String psfilestr = FileUtils.fileToString(psFile.getPath());
              TestUtils.containsAll(psfilestr,
                "Total length",
                "CNV-count",
                "CNV-percent",
                "%",
                "Sequence",
                "length",
                "CNV-count",
                "CNV-percent",
                "Lengths of CNV regions",
                "range of nt-length    : count",
                "[0               - 0] : ",
                "[1              - 10] : ",
                "[11            - 100] : ",
                "[101         - 1,000] : ",
                "[1,001      - 10,000] : ",
                "[1,0001    - 100,000] : ",
                "[..      - 1,000,000] : ",
                "[..     - 10,000,000] : ",
                "[..    - 100,000,000] : ",
                "[..  - 1,000,000,000] : ");
            }
          } finally {
            twin.close();
          }
        } finally {

          output.close();
        }
        final String cnvfilestr = FileUtils.fileToString(cnvs);
        //System.err.println(cnvfilestr);
        TestUtils.containsAll(cnvfilestr, "#Seq" + TB + "start" + TB + "end" + TB + "label" + TB + "cn" + TB + "bp-cn" + TB + "error",
          "a" + TB + "0" + TB + "35" + TB + "cnv" + TB + "2" + TB + "0" + TB + "0.0");

        final ByteArrayOutputStream bos = new ByteArrayOutputStream();

        new CnvFileChecker(new PrintStream(bos), new StringReader(cnvfilestr)).check();

        bos.close();

        assertEquals("", bos.toString());

        final int totalLength = calculateTotalLength(cnvfilestr);
        final SequencesReader outputReader = SequencesReaderFactory.createDefaultSequencesReader(outputDirectory);
        try {
          final SequencesReader twinReader = SequencesReaderFactory.createDefaultSequencesReader(twinDirectory);
          try {
            assertEquals(totalLength, outputReader.totalLength() + twinReader.totalLength());
          } finally {
            twinReader.close();
          }
        } finally {
          outputReader.close();
        }
      } finally {
        dsr.close();

      }
    } finally {
    }
  }

  public void testGenerateLengthTwoSequences() throws IOException {
    Diagnostic.setLogStream();
    //final File temp = FileHelper.createTempDirectory();
    final String genome = ""
      + ">a" + StringUtils.LS
      + "ATTGCAGCTATTGCAGCTATTGCAGCTATTGCAGCTATTGCAGCTA" + StringUtils.LS
      + ">b" + StringUtils.LS
      + "CATGATTGCAGCTAGCATCGTGTCACA" + StringUtils.LS;
    final File in = ReaderTestUtils.getDNASubDir(genome, mDir);
    try {
      final SequencesReader dsr = SequencesReaderFactory.createDefaultSequencesReader(in);
      try {
        final File cnvs = new File(mDir, "test2.cnv");
        final File outputDirectory = new File(mDir, "out2");
        final File twinDirectory = new File(mDir, "twin2");

        final SdfWriter output = new SdfWriter(outputDirectory, Constants.MAX_FILE_SIZE,
            PrereadType.UNKNOWN, false, true, false, dsr.type());
        try {
          final SdfWriter twin = new SdfWriter(twinDirectory, Constants.MAX_FILE_SIZE,
              PrereadType.UNKNOWN, false, true, false, dsr.type());
          try {
            final CnvSimulator cs = new CnvSimulator(dsr, output, twin, new FileOutputStream(cnvs), mRandom, mPriors, 10,  Integer.MAX_VALUE);
            //cs.setPriors(mPriors);
            cs.generate();
            final String csstr = cs.toString();
            //System.err.println(csstr);
            TestUtils.containsAll(csstr, "CnvSimulator", "No. breakpoints ",
              "Sequence 0", "Sequence 1");

          } finally {
            twin.close();
          }
        } finally {
          output.close();
        }
        final String cnvfilestr = FileUtils.fileToString(cnvs);
        //System.err.println(cnvfilestr);
        TestUtils.containsAll(cnvfilestr, "#Seq" + TB + "start" + TB + "end" + TB + "label" + TB + "cn" + TB + "bp-cn" + TB + "error",
          "b" + TB + "27" + TB + "27" + TB + "cnv" + TB + "0" + TB + "0" + TB + "0.0");


        final int totalLength = calculateTotalLength(cnvfilestr);
        final SequencesReader outputReader = SequencesReaderFactory.createDefaultSequencesReader(outputDirectory);
        try {
          final SequencesReader twinReader = SequencesReaderFactory.createDefaultSequencesReader(twinDirectory);
          try {
            assertEquals(totalLength, outputReader.totalLength() + twinReader.totalLength());
          } finally {
            twinReader.close();
          }
        } finally {
          outputReader.close();
        }
      } finally {
        dsr.close();

      }
    } finally {
    }

  }

  /**
   * calculate total length summing up all regions
   * @param cnvfilestr filter string
   * @return total length summing up all region lengths
   */
  public static int calculateTotalLength(String cnvfilestr) {
    final String[] lines = cnvfilestr.split(StringUtils.LS);
    int copyNumber = 0;
    int totalLength = 0;
    int lastPosition = 0;
    for (int i = 3; i < lines.length; i++) {
      final String[] parts = lines[i].split(TB);
      totalLength += (Integer.parseInt(parts[1]) - lastPosition) * copyNumber;
      lastPosition = Integer.parseInt(parts[1]);
      copyNumber = Integer.parseInt(parts[4]);
    }
    return totalLength;
  }
}
