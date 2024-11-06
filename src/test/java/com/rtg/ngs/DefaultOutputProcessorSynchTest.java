/*
 * Copyright (c) 2018. Real Time Genomics Limited.
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the
 *    distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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
