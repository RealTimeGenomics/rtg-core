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
package com.rtg.ngs.tempstage;

import java.io.ByteArrayInputStream;
import java.io.File;

import com.rtg.alignment.ActionsHelper;
import com.rtg.alignment.AlignmentResult;
import com.rtg.launcher.HashingRegion;
import com.rtg.launcher.SequenceParams;
import com.rtg.mode.DnaUtils;
import com.rtg.ngs.NgsFilterParams;
import com.rtg.ngs.NgsOutputParams;
import com.rtg.ngs.NgsOutputParamsBuilder;
import com.rtg.ngs.NgsParams;
import com.rtg.ngs.NgsParamsBuilder;
import com.rtg.ngs.ReadStatusTracker;
import com.rtg.ngs.SharedResources;
import com.rtg.ngs.blocking.ReadBlocker;
import com.rtg.pairedend.ReadStatusListener;
import com.rtg.pairedend.UnfilteredHitInfo;
import com.rtg.reader.ReaderTestUtils;
import com.rtg.util.IntegerOrPercentage;
import com.rtg.util.MaxShiftFactor;
import com.rtg.util.StringUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.MemoryPrintStream;
import com.rtg.util.test.FileHelper;

import junit.framework.TestCase;


/**
 */
public class UnfilteredTempFileWriterTest extends TestCase {

  private File mTmpDir;
  private UnfilteredTempFileWriter mUsaw;
  private MockStatusTracker mListener = null;

  @Override
  public void setUp() throws Exception {
    Diagnostic.setLogStream();
    mTmpDir = FileUtils.createTempDir("usaw", "ksdjf");
  }

  private void setUpStuff(int matedMaxScore, int unmatedMaxScore) throws Exception {

    final File template = FileUtils.createTempDir("template", "ngs", mTmpDir);
    final File left = FileUtils.createTempDir("left", "ngs", mTmpDir);
    final File right = FileUtils.createTempDir("right", "ngs", mTmpDir);

    ReaderTestUtils.getReaderDNA(TEMPLATE, template, null).close();
    ReaderTestUtils.getReaderDNA(READ_LEFT, left, null).close();
    ReaderTestUtils.getReaderDNA(READ_RIGHT, right, null).close();

    final NgsFilterParams filterparams = NgsFilterParams.builder().matedMaxMismatches(new IntegerOrPercentage(matedMaxScore)).unmatedMaxMismatches(new IntegerOrPercentage(unmatedMaxScore)).create();
    final NgsOutputParams outparams = new NgsOutputParamsBuilder().filterParams(filterparams).create();
    final NgsParams params = new NgsParamsBuilder().buildFirstParams(SequenceParams.builder().directory(left).useMemReader(true).create())
          .buildSecondParams(SequenceParams.builder().directory(right).useMemReader(true).create())
          .searchParams(SequenceParams.builder().directory(template).useMemReader(true).loadNames(true).create()).outputParams(outparams).maxFragmentLength(5)
          .gapOpenPenalty(1).gapExtendPenalty(1).substitutionPenalty(1).unknownsPenalty(0).alignerBandWidthFactor(new MaxShiftFactor(0.1))
          .create();
    mListener = new MockStatusTracker();
    mUsaw = new UnfilteredTempFileWriter(mListener, SharedResources.generateSharedResources(params), params, new ReadBlocker(55, 55), new ReadBlocker(55, 55));
  }

  @Override
  public void tearDown() {
    mUsaw = null;
    FileHelper.deleteAll(mTmpDir);
    mTmpDir = null;
    mListener = null;
  }

  static class MockStatusTracker implements ReadStatusListener {

    protected String mStatusString = "";

    @Override
    public void addStatus(int readId, int status) {
      mStatusString = mStatusString + "," + readId + ":" + ReadStatusTracker.statusToString(status);
    }
  }

  static final String TEMPLATE_STR = "acacactgcaagacaagagggcctcccacagcactctcagcccacactggtcgggggccaaagggg";
  static final String TEMPLATE = ">t" + StringUtils.LS + TEMPLATE_STR + StringUtils.LS;
  static final String LEFT_READ_AS2 =  "acacactgcaagcaagagggcctccc";
  static final String LEFT_READ_AS6 =  "acacactgcggtctgagagggcctcccac";
  static final String LEFT_READ_AS11 = "ctgcaactgttctaaagctcccacagcactct";      //starts at tpos 5
  static final String RIGHT_READ_AS0 = "cccctttggcccccgaccagtgtgggctga";
  static final String RIGHT_READ_AS8 = "cccctttttaaaaagagcagtgtgggctga";
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

  public void testUnmatedScore() throws Exception {
    final int matedMaxScore = 10;
    final int unmatedMaxScore = 5;
    setUpStuff(matedMaxScore, unmatedMaxScore);
    final UnfilteredHitInfo hit = new UnfilteredHitInfo();
    hit.setValues(true, false, 0, 0);

    assertEquals(-1, hit.score());
    mUsaw.nextTemplateId(0);
    assertTrue(mUsaw.checkUnmatedScore(hit));
    assertFalse(hit.getMatedOk());
    assertTrue(hit.getUnmatedOk());
    assertEquals(2, hit.score());
//                acacactgcaagcaagagggcctccc
//                acacactgcaagacaagagggcctcccacag
    assertEquals("============D==============", hit.alignment().getActionsString());
    assertEquals(2, hit.alignment().getScore());
    assertEquals(0, hit.alignment().getStart());
    assertEquals(",0:UNMATED_COMPUTE_ALIGNMENT_FIRST ", mListener.mStatusString);

    hit.setValues(false, true, 0, TEMPLATE_STR.length() - RIGHT_READ_AS0.length());
    assertEquals(-1, hit.score());
    assertTrue(mUsaw.checkUnmatedScore(hit));
    assertTrue(mUsaw.checkUnmatedScore(hit));   //check early termination
    assertFalse(hit.getMatedOk());
    assertTrue(hit.getUnmatedOk());
    assertEquals(0, hit.score());
//r               cccctttggcccccgaccagtgtgggctga
//t               actctcagcccacactggtcgggggccaaagggg
//tr              cccctttggcccccgaccagtgtgggctgagagt
    assertEquals("==============================", hit.alignment().getActionsString());
    assertEquals(0, hit.alignment().getScore());
    assertEquals(TEMPLATE_STR.length() - RIGHT_READ_AS0.length(), hit.alignment().getStart());
    assertEquals(",0:UNMATED_COMPUTE_ALIGNMENT_FIRST ,0:UNMATED_COMPUTE_ALIGNMENT_SECOND ", mListener.mStatusString);

    hit.setValues(true, false, 1, 0);
    assertEquals(-1, hit.score());
    assertFalse(mUsaw.checkUnmatedScore(hit));
    assertFalse(mUsaw.checkUnmatedScore(hit));  //check early termination
    assertFalse(hit.getMatedOk());
    assertFalse(hit.getUnmatedOk());
    assertTrue(hit.score() > unmatedMaxScore);
    assertNull(hit.alignment());

    hit.setValues(false, true, 4, 0);
    assertEquals(-1, hit.score());
    assertFalse(mUsaw.checkUnmatedScore(hit));
    assertFalse(hit.getMatedOk());
    assertFalse(hit.getUnmatedOk());
    assertTrue(hit.score() > unmatedMaxScore);
    assertNull(hit.alignment());
  }

  public void testCheckScores() throws Exception {
    final int matedMaxScore = 10;
    final int unmatedMaxScore = 5;
    setUpStuff(matedMaxScore, unmatedMaxScore);
    mUsaw.nextTemplateId(0);

    final UnfilteredHitInfo hit = new UnfilteredHitInfo();
    final UnfilteredHitInfo mate = new UnfilteredHitInfo();

    //r both under unmated threshold
    hit.setValues(true, false, 0, 0);
    mate.setValues(false, true, 0, TEMPLATE_STR.length() - RIGHT_READ_AS0.length());
    assertTrue(mUsaw.checkScores(hit, mate));
    assertEquals(",0:UNMATED_COMPUTE_ALIGNMENT_FIRST ,0:UNMATED_COMPUTE_ALIGNMENT_SECOND ", mListener.mStatusString);
    assertTrue(hit.getMatedOk());
    assertFalse(hit.getUnmatedOk());
    assertEquals(2, hit.score());
    assertTrue(mate.getMatedOk());
    assertFalse(mate.getUnmatedOk());
    assertEquals(0, mate.score());

    //r2 left under mated but over unmated, right under both
    hit.setValues(true, false, 1, 0);
    mate.setValues(false, true, 1, TEMPLATE_STR.length() - RIGHT_READ_AS0.length());
    assertTrue(mUsaw.checkScores(hit, mate));
    assertEquals(",0:UNMATED_COMPUTE_ALIGNMENT_FIRST ,0:UNMATED_COMPUTE_ALIGNMENT_SECOND ,1:UNMATED_COMPUTE_ALIGNMENT_FIRST ,1:UNMATED_COMPUTE_ALIGNMENT_SECOND ", mListener.mStatusString);
    assertTrue(mUsaw.checkScores(hit, mate));
    assertNotNull(hit.alignment());
    assertNotNull(mate.alignment());
    assertTrue(hit.getMatedOk());
    assertFalse(hit.getUnmatedOk());
    assertEquals(6, hit.score());
    assertTrue(mate.getMatedOk());
    assertFalse(mate.getUnmatedOk());
    assertEquals(0, mate.score());

    //r3 both under mated but over unmated
    hit.setValues(true, false, 2, 0);
    mate.setValues(false, true, 2, TEMPLATE_STR.length() - RIGHT_READ_AS8.length());
    assertTrue(mUsaw.checkScores(hit, mate));
    assertTrue(hit.getMatedOk());
    assertFalse(hit.getUnmatedOk());
    assertEquals(6, hit.score());
    assertTrue(mate.getMatedOk());
    assertFalse(mate.getUnmatedOk());
    assertEquals(8, mate.score());

    //r4 left under unmated, right under mated but over unmated
    hit.setValues(true, false, 3, 0);
    mate.setValues(false, true, 3, TEMPLATE_STR.length() - RIGHT_READ_AS8.length());
    assertTrue(mUsaw.checkScores(hit, mate));
    assertTrue(hit.getMatedOk());
    assertFalse(hit.getUnmatedOk());
    assertEquals(2, hit.score());
    assertTrue(mate.getMatedOk());
    assertFalse(mate.getUnmatedOk());
    assertEquals(8, mate.score());

    //r5 left under mated but over unmated, right over both
    hit.setValues(true, false, 4, 0);
    mate.setValues(false, true, 4, 5);
    assertFalse(mUsaw.checkScores(hit, mate));
    assertFalse(hit.getMatedOk());
    assertFalse(hit.getUnmatedOk());
    assertEquals(6, hit.score());
    assertFalse(mate.getMatedOk());
    assertFalse(mate.getUnmatedOk());
    assertTrue(mate.score() > matedMaxScore);

    //r6 left over both, right under mated but over unmated
    hit.setValues(true, false, 5, 0);
    mate.setValues(false, true, 5, TEMPLATE_STR.length() - RIGHT_READ_AS8.length());
    assertFalse(mUsaw.checkScores(hit, mate));
    assertFalse(hit.getMatedOk());
    assertFalse(hit.getUnmatedOk());
    assertTrue(hit.score() > matedMaxScore);
    assertFalse(mate.getMatedOk());
    assertFalse(mate.getUnmatedOk());
    assertEquals(-1, mate.score());     //we want to die early if the left hit exceeds both thresholds

    //r7 left under both, right over both
    hit.setValues(true, false, 6, 0);
    mate.setValues(false, true, 6, 5);
    assertTrue(mUsaw.checkScores(hit, mate));
    assertFalse(hit.getMatedOk());
    assertTrue(hit.getUnmatedOk());
    assertEquals(2, hit.score());
    assertFalse(mate.getMatedOk());
    assertFalse(mate.getUnmatedOk());
    assertTrue(mate.score() > matedMaxScore);

    //r8 left over both, right under both
    hit.setValues(true, false, 7, 0);
    mate.setValues(false, true, 7, TEMPLATE_STR.length() - RIGHT_READ_AS0.length());
    assertFalse(mUsaw.checkScores(hit, mate));
    assertFalse(hit.getMatedOk());
    assertFalse(hit.getUnmatedOk());
    assertTrue(hit.score() > matedMaxScore);
    assertFalse(mate.getMatedOk());
    assertFalse(mate.getUnmatedOk());
    assertEquals(-1, mate.score());     //we want to die early if the left hit exceeds both thresholds
  }

  public void testCheckScore2() throws Exception {
    final int matedMaxScore = 5;
    final int unmatedMaxScore = 10;
    setUpStuff(matedMaxScore, unmatedMaxScore);
    mUsaw.nextTemplateId(0);

    final UnfilteredHitInfo hit = new UnfilteredHitInfo();
    final UnfilteredHitInfo mate = new UnfilteredHitInfo();

    //r2 left under unmated but over mated, right under both but wont be calculated yet
    hit.setValues(true, false, 1, 0);
    mate.setValues(false, true, 1, TEMPLATE_STR.length() - RIGHT_READ_AS0.length());
    assertTrue(mUsaw.checkScores(hit, mate));
    assertFalse(hit.getMatedOk());
    assertTrue(hit.getUnmatedOk());
    assertEquals(6, hit.score());
    assertFalse(mate.getMatedOk());
    assertFalse(mate.getUnmatedOk());
    assertEquals(-1, mate.score());

    //r1 both under, right first
    hit.setValues(false, true, 0, TEMPLATE_STR.length() - RIGHT_READ_AS0.length());
    mate.setValues(true, false, 0, 0);
    assertTrue(mUsaw.checkScores(hit, mate));
    assertTrue(hit.getMatedOk());
    assertFalse(hit.getUnmatedOk());
    assertEquals(0, hit.score());
    assertTrue(mate.getMatedOk());
    assertFalse(mate.getUnmatedOk());
    assertEquals(2, mate.score());
    assertEquals(",1:UNMATED_COMPUTE_ALIGNMENT_FIRST ,0:UNMATED_COMPUTE_ALIGNMENT_SECOND ,0:UNMATED_COMPUTE_ALIGNMENT_FIRST ", mListener.mStatusString);

    //r4 left under mated, right under unmated but over mated
    hit.setValues(true, false, 3, 0);
    mate.setValues(false, true, 3, TEMPLATE_STR.length() - RIGHT_READ_AS8.length());
    assertTrue(mUsaw.checkScores(hit, mate));
    assertFalse(hit.getMatedOk());
    assertTrue(hit.getUnmatedOk());
    assertEquals(2, hit.score());
    assertFalse(mate.getMatedOk());
    assertTrue(mate.getUnmatedOk());
    assertEquals(8, mate.score());
  }

  public void testWriting() throws Exception {
    setUpStuff(10, 5);
    final MemoryPrintStream mps = new MemoryPrintStream();

    mUsaw.initialiseUnmated(mps.outputStream());

    final UnfilteredHitInfo h = new UnfilteredHitInfo();
    h.setValues(true, false, 0, 0);
    h.setAlignment(2, new AlignmentResult(DnaUtils.encodeString(LEFT_READ_AS2), ActionsHelper.build("============D==============", 0, 2), DnaUtils.encodeString(TEMPLATE_STR)));
    h.setUnmatedOk(true);

    mUsaw.nextTemplateId(0);
    mUsaw.unmatedResultUnfiltered(h);
    assertFalse(h.alignment().isReverse());
    h.setUnmatedOk(false);

    h.setValues(false, true, 0, TEMPLATE_STR.length() - RIGHT_READ_AS8.length());
    h.setAlignment(8, new AlignmentResult(DnaUtils.encodeString(RIGHT_READ_AS8), ActionsHelper.build("==============================", TEMPLATE_STR.length() - RIGHT_READ_AS8.length(), 8), DnaUtils.encodeString(TEMPLATE_STR)));
    h.setMatedOk(true);
    mUsaw.unmatedResultUnfiltered(h);
    assertTrue(h.alignment().isReverse());

    mUsaw.setClipRegion(new HashingRegion(99999, 999999));
    mUsaw.unmatedResultUnfiltered(h);

    mUsaw.close();
    final TempRecordReaderNio dis = new TempRecordReaderNio(new ByteArrayInputStream(mps.toByteArray()), new TempRecordReader.RecordFactory(true, false, false, true));
    BinaryTempFileRecord rec = dis.readRecord();
    assertNotNull(rec);
    assertEquals(0, rec.getReadId());
    assertEquals(65, rec.getSamFlags());
    assertEquals(0, rec.getReferenceId());
    assertEquals(1, rec.getStartPosition());
    assertEquals("12=1D14=", new String(rec.getCigarString()));
    assertEquals(2, rec.getAlignmentScore());
    assertEquals(1, rec.getNumberMismatches());

    rec = dis.readRecord();
    assertNotNull(rec);
    assertEquals(0, rec.getReadId());
    assertEquals(145, rec.getSamFlags() & 0xff);
    assertEquals(0, rec.getReferenceId());
    assertEquals(37, rec.getStartPosition());
    assertEquals("30=", new String(rec.getCigarString()));
    assertEquals(8, rec.getAlignmentScore());
    assertEquals(0, rec.getNumberMismatches());

    assertEquals(",0:UNMATED_ALIGN_SCORE_FIRST ,0:UNMATED_FIRST ,0:UNMATED_ALIGN_SCORE_SECOND ,0:UNMATED_SECOND ", mListener.mStatusString);

  }
}
