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

package com.rtg.ngs;

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.IOException;

import com.rtg.launcher.SequenceParams;
import com.rtg.reader.ReaderTestUtils;
import com.rtg.util.IntegerOrPercentage;
import com.rtg.util.StringUtils;
import com.rtg.util.cli.CommandLine;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.MemoryPrintStream;
import com.rtg.util.test.FileHelper;

import htsjdk.samtools.SAMRecord;

import junit.framework.TestCase;

/**
 */
public class UnmappedSamAlignmentWriterTest extends TestCase {

  private File mDir;

  @Override
  public void setUp() throws Exception {
    mDir = FileHelper.createTempDirectory();
    Diagnostic.setLogStream();
  }

  @Override
  public void tearDown() {
    FileHelper.deleteAll(mDir);
    mDir = null;
  }


  private static NgsParams getCommonTestParams(File left, File right, File template, IntegerOrPercentage maxMatedScore, IntegerOrPercentage maxUnmatedScore) {
    final NgsOutputParams op = NgsOutputParams.builder().filterParams(NgsFilterParams.builder()
        .matedMaxMismatches(maxMatedScore)
        .unmatedMaxMismatches(maxUnmatedScore)
        .outputFilter(OutputFilter.PAIRED_END).maxTopResults(10).create()).create();
    return NgsParams.builder().buildFirstParams(SequenceParams.builder().directory(left).useMemReader(true).create()).buildSecondParams(SequenceParams.builder().directory(right).useMemReader(true).create()).searchParams(SequenceParams.builder().directory(template).loadNames(true).useMemReader(true).create()).maxFragmentLength(1000).minFragmentLength(0).outputParams(op).create();
  }

  static final String TEMPLATE = ">t" + StringUtils.LS + "acacactgcaagacaagagggcctcccacagcactctcagcccacactggtcgggggccaaagggg" + StringUtils.LS;
  static final String TEMP_LEFT = "acacactgcaagcaagagggcctccc";
  static final String TEMP_RIGHT = "cccctttggcccccgaccagtgtgggctga";
  static final String READ_LEFT = ">r" + StringUtils.LS + TEMP_LEFT + StringUtils.LS;
  static final String READ_RIGHT = ">r" + StringUtils.LS + TEMP_RIGHT + StringUtils.LS;

  public void testPairing() throws IOException {
    final File template = FileUtils.createTempDir("template", "ngs", mDir);
    final File left = FileUtils.createTempDir("left", "ngs", mDir);
    final File right = FileUtils.createTempDir("right", "ngs", mDir);

    ReaderTestUtils.getReaderDNA(TEMPLATE, template, null).close();
    ReaderTestUtils.getReaderDNA(READ_LEFT, left, null).close();
    ReaderTestUtils.getReaderDNA(READ_RIGHT, right, null).close();

    CommandLine.setCommandArgs("wibble", "-h", "dream-turnip");
    try {
      final NgsOutputParams op = NgsOutputParams.builder().filterParams(NgsFilterParams.builder().outputFilter(OutputFilter.PAIRED_END).create()).create();
      final NgsParams param = NgsParams.builder().buildFirstParams(SequenceParams.builder().directory(left).useMemReader(true).create()).buildSecondParams(SequenceParams.builder().directory(right).useMemReader(true).create()).searchParams(SequenceParams.builder().directory(template).loadNames(true).useMemReader(true).create()).maxFragmentLength(1000).minFragmentLength(0).outputParams(op).create();

      final SharedResources sr = SharedResources.generateSharedResources(param);
      final ByteArrayOutputStream bos = new ByteArrayOutputStream();
      try {
        final UnmappedSamAlignmentWriter w = new UnmappedSamAlignmentWriter(param.outputParams().tempFilesDirectory(), sr.getHeader());
        UnmappedSamRecordFactory fact = new UnmappedSamRecordFactory(param, sr);
        try {
          w.initialiseUnmapped(bos, false, true, true);
          w.unmappedRecord(fact.unmappedResult(0, true, (char) 1, false));
        } finally {
          w.close();
        }
      } finally {
        bos.close();
      }
      final String us = bos.toString();

      assertTrue(us.contains("@PG"));
      assertTrue(us.contains("CL:wibble -h dream-turnip"));
      assertTrue(us.contains("@SQ"));
      assertTrue(us.contains("0\t69\t*\t0\t0\t*\t*\t0\t0\tACACACTGCAAGCAAGAGGGCCTCCC\t*\tXC:A:"));

    } finally {
      CommandLine.clearCommandArgs();
    }
  }


  static final String TEMPLATE_3 = ">t" + StringUtils.LS + "tttaccccccccccccc" + StringUtils.LS;
  static final String TEMP_LEFT_3 = "ttt";
  static final String TEMP_RIGHT_3 = "ccc";
  static final String READ_LEFT_3 = ">r" + StringUtils.LS + TEMP_LEFT_3 + StringUtils.LS;
  static final String READ_RIGHT_3 = ">r" + StringUtils.LS + TEMP_RIGHT_3 + StringUtils.LS;

  public void testOtherStuff() throws Exception {
    final File template = FileUtils.createTempDir("template", "ngs", mDir);
    final File left = FileUtils.createTempDir("left", "ngs", mDir);
    final File right = FileUtils.createTempDir("right", "ngs", mDir);

    ReaderTestUtils.getReaderDNA(READ_LEFT_3, left, null).close();
    ReaderTestUtils.getReaderDNA(READ_RIGHT_3, right, null).close();

    final String templ = TEMPLATE_3;

    ReaderTestUtils.getReaderDNA(templ, template, null).close();
    final NgsParams param = getCommonTestParams(left, right, template, IntegerOrPercentage.valueOf(5), IntegerOrPercentage.valueOf(5));

    final MemoryPrintStream baos = new MemoryPrintStream();
    final ByteArrayOutputStream baosunmapped = new ByteArrayOutputStream();
    Diagnostic.setLogStream(baos.printStream());
    try {
      final SharedResources sr = SharedResources.generateSharedResources(param);
      final UnmappedSamAlignmentWriter w = new UnmappedSamAlignmentWriter(param.outputParams().tempFilesDirectory(), sr.getHeader());
      final UnmappedSamRecordFactory fact = new UnmappedSamRecordFactory(param, sr);
      try {
        final SAMRecord rec = fact.createUnmappedSamRecord(0, "", new byte[0], true, false, 'x', false);
        assertTrue(rec.getSecondOfPairFlag());
        w.initialiseUnmapped(baosunmapped, false, true, true);
        try {
          w.unmappedRecord(fact.unmappedResult(0, true, '\0', false));
        } finally {
          w.close();
        }
      } finally {
        fact.close();
      }
      assertNull(fact.mReader1);
      assertNull(fact.mReader2);

      final String unmappedstring = baosunmapped.toString();
      assertTrue(unmappedstring.contains("0\t69\t*"));
      assertFalse(unmappedstring.contains("0\t67\t0"));
      assertFalse(unmappedstring.contains("0\t89\t0"));
      assertFalse(unmappedstring.contains("0\t131\t0"));
    } finally {
      Diagnostic.setLogStream();
    }
  }
}
