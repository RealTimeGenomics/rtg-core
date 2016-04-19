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
package com.rtg.pairedend;

import java.io.File;
import java.io.OutputStream;

import com.rtg.launcher.SequenceParams;
import com.rtg.ngs.NgsFilterParams;
import com.rtg.ngs.NgsOutputParams;
import com.rtg.ngs.NgsOutputParamsBuilder;
import com.rtg.ngs.NgsParams;
import com.rtg.ngs.NgsParamsBuilder;
import com.rtg.ngs.SharedResources;
import com.rtg.ngs.blocking.ReadBlocker;
import com.rtg.ngs.tempstage.PairedTempFileWriterImplTest;
import com.rtg.ngs.tempstage.UnfilteredTempFileWriter;
import com.rtg.reader.ReaderTestUtils;
import com.rtg.util.IntegerOrPercentage;
import com.rtg.util.StringUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.io.FileUtils;
import com.rtg.util.machine.MachineOrientation;
import com.rtg.util.test.FileHelper;

import junit.framework.TestCase;

/**
 */
public class UnfilteredSlidingWindowCollectorTest extends TestCase {

  private class TestUnfilteredAlignmentWriter extends UnfilteredTempFileWriter {

    TestUnfilteredAlignmentWriter(ReadStatusListener listener, SharedResources sr, NgsParams params) {
      super(listener, sr, params, new ReadBlocker(55, 55), new ReadBlocker(55, 55));
      //mExpectedTemplateCount = expectedTemplateCount;
      //mExpectedPairResultCount = expectedPairResultCount;
    }

    @Override
    public void initialiseUnmated(OutputStream unmatedOut) {
    }
    @Override
    public boolean checkUnmatedScore(UnfilteredHitInfo hit) {
      if (hit.equals(mUnhit)) {
        return false;
      }
      assertEquals(mMate, hit);
      mMate.setAlignment(0, null);
      mMate.setUnmatedOk(true);
      return true;
    }

    @Override
    public boolean checkScores(UnfilteredHitInfo hit, UnfilteredHitInfo mate) {
      if (hit.equals(mUnhit)) {
        return false;
      }
      assertEquals(mHit, hit);
      mHit.setAlignment(2, null);
      mHit.setMatedOk(true);
      assertEquals(mMate, mate);
      return true;
    }

    @Override
    public void unmatedResultUnfiltered(UnfilteredHitInfo hit) {
      if (hit.getMatedOk()) {
        assertEquals(mHit, hit);
      } else {
        assertEquals(mMate, hit);
      }
    }
  }


  private File mTmpDir = null;
  private NgsParams mParams = null;
  private static final int UNMATED_MAX_SCORE = 5;
  private static final int MATED_MAX_SCORE = 10;

  /*
   * r = as 2, 0
   * r2 = as 6, 0
   * r3 = as 6, 8
   * r4 = as 2, 8
   * r5 = as 6, 11
   * r6 = as 11, 8
   * r7 = as 2, 11
   * r8 = as 11, 0
   */
  static final String TEMPLATE_STR = "acacactgcaagacaagagggcctcccacagcactctcagcccacactggtcgggggccaaagggg";
  static final String TEMPLATE = ">t" + StringUtils.LS + TEMPLATE_STR + StringUtils.LS;
  static final String LEFT_READ_AS2 =  "acacactgcaagcaagagggcctccc";
  static final String LEFT_READ_AS6 =  "acacactgcggtctgagagggcctcccac";
  static final String RIGHT_READ_AS0 = "cccctttggcccccgaccagtgtgggctga";
  static final String RIGHT_READ_AS8 = "cccctttttaaaaagagcagtgtgggctga";
  static final String LEFT_READ_AS11 = "ctgcaactgttctaaagctcccacagcactct";      //starts at tpos 5
  static final String RIGHT_READ_AS11 = "agagtgctgtgggagctttagaacagttgcag";      //starts at tpos 5
  static final String READ_LEFT = ">r" + StringUtils.LS + LEFT_READ_AS2 + StringUtils.LS
                          + ">r2" + StringUtils.LS + LEFT_READ_AS6 + StringUtils.LS
                          + ">r3" + StringUtils.LS + LEFT_READ_AS6 + StringUtils.LS
                          + ">r4" + StringUtils.LS + LEFT_READ_AS2 + StringUtils.LS
                          + ">r5" + StringUtils.LS + LEFT_READ_AS6 + StringUtils.LS
                          + ">r6" + StringUtils.LS + LEFT_READ_AS11 + StringUtils.LS
                          + ">r7" + StringUtils.LS + LEFT_READ_AS2 + StringUtils.LS
                          + ">r8" + StringUtils.LS + LEFT_READ_AS11 + StringUtils.LS;
  static final String READ_RIGHT = ">r" + StringUtils.LS + RIGHT_READ_AS0 + StringUtils.LS
                          + ">r2" + StringUtils.LS + RIGHT_READ_AS0 + StringUtils.LS
                          + ">r3" + StringUtils.LS + RIGHT_READ_AS8 + StringUtils.LS
                          + ">r4" + StringUtils.LS + RIGHT_READ_AS8 + StringUtils.LS
                          + ">r5" + StringUtils.LS + RIGHT_READ_AS11 + StringUtils.LS
                          + ">r6" + StringUtils.LS + RIGHT_READ_AS8 + StringUtils.LS
                          + ">r7" + StringUtils.LS + RIGHT_READ_AS11 + StringUtils.LS
                          + ">r8" + StringUtils.LS + RIGHT_READ_AS0 + StringUtils.LS;

  @Override
  public void setUp() throws Exception {
    Diagnostic.setLogStream();
    final File tmpDir = FileUtils.createTempDir("usaw", "ksdjf");

    final File template = FileUtils.createTempDir("template", "ngs", tmpDir);
    final File left = FileUtils.createTempDir("left", "ngs", tmpDir);
    final File right = FileUtils.createTempDir("right", "ngs", tmpDir);
    //final File out = File.createTempFile("sam", "out", tmpDir);

    ReaderTestUtils.getReaderDNA(TEMPLATE, template, null).close();
    ReaderTestUtils.getReaderDNA(READ_LEFT, left, null).close();
    ReaderTestUtils.getReaderDNA(READ_RIGHT, right, null).close();

    final NgsFilterParams filterparams = NgsFilterParams.builder().matedMaxMismatches(new IntegerOrPercentage(MATED_MAX_SCORE)).unmatedMaxMismatches(new IntegerOrPercentage(UNMATED_MAX_SCORE)).create();
    final NgsOutputParams outparams = new NgsOutputParamsBuilder().filterParams(filterparams).create();
    mParams = new NgsParamsBuilder().buildFirstParams(SequenceParams.builder().directory(left).useMemReader(true).create())
          .buildSecondParams(SequenceParams.builder().directory(right).useMemReader(true).create())
          .searchParams(SequenceParams.builder().directory(template).useMemReader(true).loadNames(true).create()).outputParams(outparams).maxFragmentLength(5).create();
    mTmpDir = tmpDir;
  }

  @Override
  public void tearDown() {
    mParams = null;
    mHit = null;
    mMate = null;
    mUnhit = null;
    FileHelper.deleteAll(mTmpDir);
    mTmpDir = null;
  }

  UnfilteredHitInfo mHit = null;
  UnfilteredHitInfo mMate = null;
  UnfilteredHitInfo mUnhit = null;

  public void testFlush() throws Exception {
    final SharedResources sr = SharedResources.generateSharedResources(mParams);
    final UnfilteredTempFileWriter mUsaw = new TestUnfilteredAlignmentWriter(new PairedTempFileWriterImplTest.UselessStatusIdListener(), sr, mParams);
    final UnfilteredSlidingWindowCollector uswc = new UnfilteredSlidingWindowCollector(70, 0, MachineOrientation.ANY, mUsaw, sr, mParams.outputParams().calibrateRegions());

    uswc.nextTemplateId(0);
    uswc.mCurrentReferencePosition = -1;

    mHit = uswc.createHitInfo();
    mMate = uswc.createHitInfo();
    mUnhit = uswc.createHitInfo();
    assertNotNull(mHit);
    assertTrue(uswc.checkPair(mHit, mMate));

    mHit.setValues(true, false, 0, 0);
    uswc.setReadLookup(uswc.readHash(0, mHit.first()), mHit);
    uswc.mReadsWindow[0].add(mHit);
    uswc.mReadsWindowInUse[0]++;

    mMate.setValues(false, true, 0, TEMPLATE_STR.length() - RIGHT_READ_AS0.length());
    uswc.mReadsWindow[TEMPLATE_STR.length() - RIGHT_READ_AS0.length()].add(mMate);
    uswc.mReadsWindowInUse[TEMPLATE_STR.length() - RIGHT_READ_AS0.length()]++;
    uswc.setReadLookup(uswc.readHash(0, !mHit.first()), mMate);

    mUnhit.setValues(true, false, 0, 0);
    uswc.mReadsWindow[0].add(mUnhit);
    uswc.mReadsWindowInUse[0]++;
    mHit.setNext(mUnhit);

    uswc.flushToPosition(0);
    assertEquals(0, uswc.mReadsWindowInUse[0]);

    uswc.flushToPosition(TEMPLATE_STR.length());

    assertTrue(mHit.getMatedOk());
    assertFalse(mHit.getUnmatedOk());
    assertFalse(mMate.getMatedOk());
    assertTrue(mMate.getUnmatedOk());
    assertFalse(mUnhit.getMatedOk());
    assertFalse(mUnhit.getUnmatedOk());
  }

}
