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


import static com.rtg.util.StringUtils.LS;

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;

import com.rtg.index.hash.ngs.OutputProcessor;
import com.rtg.launcher.HashingRegion;
import com.rtg.util.io.FileUtils;
import com.rtg.util.test.FileHelper;

import junit.framework.TestCase;

/**
 */
public class DefaultOutputProcessorSynchTest extends TestCase {

  public void testProcessor() throws IOException {
    final File dir = FileUtils.createTempDir("dop", "test");
    try {
      final ByteArrayOutputStream b = new ByteArrayOutputStream();
      try {
        final DefaultOutputProcessorSynch dop = new DefaultOutputProcessorSynch(NgsParams.builder().outputParams(new NgsTestUtils.OverriddenNgsOutputParams(NgsTestUtils.OverriddenNgsOutputParams.builder().outStream(b).outputDir(dir))).create());
        final OutputProcessor op = dop.threadClone(HashingRegion.NONE);
        op.process(1234, "R", 123456, 30000, 2, 0);
        op.threadFinish();
        dop.finish();
      } finally {
        b.close();
      }
      assertEquals("#template-id\tframe\tread-id\ttemplate-start\tscore\tscore-indel" + LS + "1234\tR\t123456\t30001\t2\t0" + LS, b.toString());
    } finally {
      assertTrue(FileHelper.deleteAll(dir));
    }
  }

  public void testThreadClone() throws IOException {
    final File dir = FileUtils.createTempDir("dop", "test");
    try {
      final ByteArrayOutputStream b = new ByteArrayOutputStream();
      final DefaultOutputProcessorSynch dop = new DefaultOutputProcessorSynch(NgsParams.builder().outputParams(new NgsTestUtils.OverriddenNgsOutputParams(NgsTestUtils.OverriddenNgsOutputParams.builder().outStream(b).outputDir(dir))).create());
      assertTrue(dop.threadClone(HashingRegion.NONE) != dop);
      dop.close();
    } finally {
      assertTrue(FileHelper.deleteAll(dir));
    }
  }

  private static final String OUT_NOTHING = "";
  private static final String SOMETHING = "edfg\n";
  private static final String OUT_SOMETHING = "#abcd\n" + SOMETHING;

  public void testAppend() throws IOException {
    final File temp = FileUtils.createTempDir("dop", "test");
    try {
      ByteArrayOutputStream out = new ByteArrayOutputStream();
      final ArrayList<File> tempFiles = new ArrayList<>();
      tempFiles.add(new File(temp, DefaultOutputProcessorSynch.THREAD_FILE_NAME + 0));
      tempFiles.add(new File(temp, DefaultOutputProcessorSynch.THREAD_FILE_NAME + 1));
      tempFiles.add(new File(temp, DefaultOutputProcessorSynch.THREAD_FILE_NAME + 2));
      tempFiles.add(new File(temp, DefaultOutputProcessorSynch.THREAD_FILE_NAME + 3));

      FileUtils.stringToFile(OUT_NOTHING, tempFiles.get(0));
      FileUtils.stringToFile(OUT_NOTHING, tempFiles.get(1));
      FileUtils.stringToFile(OUT_NOTHING, tempFiles.get(2));
      FileUtils.stringToFile(OUT_SOMETHING, tempFiles.get(3));
      DefaultOutputProcessorSynch.appendAllFiles(out, tempFiles.toArray(new File[4]), false);
      assertEquals(OUT_SOMETHING, out.toString());
      out = new ByteArrayOutputStream();
      tempFiles.add(new File(temp, DefaultOutputProcessorSynch.THREAD_FILE_NAME + 4));
      FileUtils.stringToFile(OUT_SOMETHING, tempFiles.get(4));
      DefaultOutputProcessorSynch.appendAllFiles(out, tempFiles.toArray(new File[5]), false);
      assertEquals(OUT_SOMETHING + SOMETHING, out.toString());
      tempFiles.add(new File(temp, DefaultOutputProcessorSynch.THREAD_FILE_NAME + 5));
      tempFiles.add(new File(temp, DefaultOutputProcessorSynch.THREAD_FILE_NAME + 6));
      FileUtils.stringToFile(OUT_NOTHING, tempFiles.get(5));
      FileUtils.stringToFile(OUT_SOMETHING, tempFiles.get(6));
      out = new ByteArrayOutputStream();
      DefaultOutputProcessorSynch.appendAllFiles(out, tempFiles.toArray(new File[7]), false);
      assertEquals(OUT_SOMETHING + SOMETHING + SOMETHING, out.toString());
    } finally {
      assertTrue(FileHelper.deleteAll(temp));
    }
  }

  public void testSortedAppend() throws IOException {
    final File dir = FileUtils.createTempDir("dop", "test");
    try {
      final ByteArrayOutputStream b = new ByteArrayOutputStream();
      try {
        final DefaultOutputProcessorSynch dop = new DefaultOutputProcessorSynch(NgsParams.builder().outputParams(new NgsTestUtils.OverriddenNgsOutputParams(NgsTestUtils.OverriddenNgsOutputParams.builder().outStream(b).outputDir(dir))).create());
        final OutputProcessor one = dop.threadClone(new HashingRegion(2, 0, 2, 10, -1, -1));
        final OutputProcessor two = dop.threadClone(new HashingRegion(3, 0, 3, 10, -1, -1));
        final OutputProcessor three = dop.threadClone(new HashingRegion(1, 0, 1, 10, -1, -1));
        one.process(2, "F", 2, 2, 0, 0);
        two.process(3, "F", 3, 2, 0, 0);
        three.process(1, "F", 4, 2, 0, 0);
        one.threadFinish();
        two.threadFinish();
        three.threadFinish();
        dop.finish();
      } finally {
        b.close();
      }
      assertEquals("#template-id\tframe\tread-id\ttemplate-start\tscore\tscore-indel" + LS
              + "1\tF\t4\t3\t0\t0" + LS
              + "2\tF\t2\t3\t0\t0" + LS
              + "3\tF\t3\t3\t0\t0" + LS
              , b.toString());
    } finally {
      assertTrue(FileHelper.deleteAll(dir));
    }

  }
}
