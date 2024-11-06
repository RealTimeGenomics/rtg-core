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
package com.rtg.variant.cnv;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.Map;

import com.rtg.launcher.OutputParams;
import com.rtg.sam.DefaultSamFilter;
import com.rtg.sam.RecordIterator;
import com.rtg.sam.SamFilterIterator;
import com.rtg.sam.ThreadedMultifileIteratorTest;
import com.rtg.util.StringUtils;
import com.rtg.util.TestUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.MemoryPrintStream;
import com.rtg.util.io.TestDirectory;
import com.rtg.util.test.FileHelper;

import htsjdk.samtools.SAMRecord;

import junit.framework.TestCase;

/**
 *
 */
public class CnvProductTaskTest extends TestCase {

  @Override
  public void setUp() {
    Diagnostic.setLogStream();
  }

  public void testSequenceLengths() throws IOException {
    final File dir = FileUtils.createTempDir("cnvproducttask2", "test");
    try {
      try (RecordIterator<SAMRecord> multiIt = ThreadedMultifileIteratorTest.getIterator(dir, "com/rtg/variant/cnv/resources/testFilter.sam")) {
        final int[] lengths = CnvProductTask.makeTemplateLengths(multiIt.header().getSequenceDictionary());
        assertEquals(3, lengths.length);
        assertEquals(143000, lengths[0]);
        assertEquals(154000, lengths[1]);
        assertEquals(100000, lengths[2]);
      }
    } finally {
      assertTrue(FileHelper.deleteAll(dir));
    }
  }
  public void testSequenceNames() throws IOException {
    final File dir = FileUtils.createTempDir("cnvproducttask3", "test");
    try {
      try (RecordIterator<SAMRecord> multiIt = ThreadedMultifileIteratorTest.getIterator(dir, "com/rtg/variant/cnv/resources/testFilter.sam")) {
        final Map<Integer, String> map = CnvProductTask.makeTemplateNameMap(multiIt.header().getSequenceDictionary());
        assertEquals(3, map.size());
        assertEquals("simulatedSequence1", map.get(0));
        assertEquals("simulatedSequence2", map.get(1));
        assertEquals("simulatedSequence3", map.get(2));
      }
    } finally {
      assertTrue(FileHelper.deleteAll(dir));
    }
  }

  public void testFillBuckets() {
    final int[] startPosCounts = {0, 2, 1, 1, 0, 0, 0, 2, 1, 1, 0};
    CnvProductParams params = CnvProductParams.builder().filterStartPositions(false).bucketSize(3).create();
    CnvProductTask task = new CnvProductTask(params, TestUtils.getNullOutputStream());
    int[][] buckets = {new int[] {0, 0, 0, 0}};
    task.fillBucket(startPosCounts, buckets, 0);
    assertEquals(3, buckets[0][0]);
    assertEquals(1, buckets[0][1]);
    assertEquals(3, buckets[0][2]);
    assertEquals(1, buckets[0][3]);

    params = CnvProductParams.builder().filterStartPositions(true).bucketSize(3).create();
    task = new CnvProductTask(params, TestUtils.getNullOutputStream());
    buckets = new int[][] {new int[] {0, 0, 0, 0}};
    task.fillBucket(startPosCounts, buckets, 0);
    assertEquals(2, buckets[0][0]);
    assertEquals(1, buckets[0][1]);
    assertEquals(2, buckets[0][2]);
    assertEquals(1, buckets[0][3]);

    buckets = new int[][] {new int[] {0, 0, 0, 0}};
    task.fillBucket(startPosCounts, buckets, -1);
    assertEquals(0, buckets[0][0]);
    assertEquals(0, buckets[0][1]);
    assertEquals(0, buckets[0][2]);
    assertEquals(0, buckets[0][3]);

    buckets = new int[][] {new int[] {0, 0, 0, 0}};
    task.fillBucket(null, buckets, 0);
    assertEquals(0, buckets[0][0]);
    assertEquals(0, buckets[0][1]);
    assertEquals(0, buckets[0][2]);
    assertEquals(0, buckets[0][3]);
  }

  public void testChunkSamFile() throws IOException {
    final File dir = FileUtils.createTempDir("cnvproducttask4", "test");
    try {
      final CnvProductParams params = CnvProductParams.builder().bucketSize(6).create();
      try (SamFilterIterator it = new SamFilterIterator(ThreadedMultifileIteratorTest.getIterator(dir, "com/rtg/variant/cnv/resources/babySam.sam"), new DefaultSamFilter(params.filterParams()))) {
        final CnvProductTask task = new CnvProductTask(params, TestUtils.getNullOutputStream());
        task.setSequenceLengths(CnvProductTask.makeTemplateLengths(it.header().getSequenceDictionary()));
        final int[][] buckets = task.chunkSamFile(it);
        assertEquals(3, buckets.length);
        assertEquals(5, buckets[0].length);
        assertEquals(5, buckets[1].length);
        assertEquals(5, buckets[2].length);
        assertEquals(3, buckets[0][0]);
        assertEquals(2, buckets[0][1]);
        assertEquals(2, buckets[0][2]);
        assertEquals(0, buckets[0][3]);
        assertEquals(0, buckets[0][4]);
        assertEquals(4, buckets[1][0]);
        assertEquals(3, buckets[1][1]);
        assertEquals(1, buckets[1][2]);
        assertEquals(0, buckets[1][3]);
        assertEquals(0, buckets[1][4]);
        assertEquals(3, buckets[2][0]);
        assertEquals(3, buckets[2][1]);
        assertEquals(1, buckets[2][2]);
        assertEquals(0, buckets[2][3]);
        assertEquals(0, buckets[2][4]);
      }
    } finally {
      assertTrue(FileHelper.deleteAll(dir));
    }
  }
  public void testChunkSamFile2() throws IOException {
    final File dir = FileUtils.createTempDir("cnvproducttask5", "test");
    try {
      final CnvProductParams params = CnvProductParams.builder().bucketSize(5).create();
      try (SamFilterIterator it = new SamFilterIterator(ThreadedMultifileIteratorTest.getIterator(dir, "com/rtg/variant/cnv/resources/babySam.sam"), new DefaultSamFilter(params.filterParams()))) {
        final CnvProductTask task = new CnvProductTask(params, TestUtils.getNullOutputStream());
        task.setSequenceLengths(CnvProductTask.makeTemplateLengths(it.header().getSequenceDictionary()));
        final int[][] buckets = task.chunkSamFile(it);
        assertEquals(3, buckets.length);
        assertEquals(5, buckets[0].length);
        assertEquals(5, buckets[1].length);
        assertEquals(5, buckets[2].length);
        assertEquals(3, buckets[0][0]);
        assertEquals(1, buckets[0][1]);
        assertEquals(2, buckets[0][2]);
        assertEquals(1, buckets[0][3]);
        assertEquals(0, buckets[0][4]);
        assertEquals(3, buckets[1][0]);
        assertEquals(3, buckets[1][1]);
        assertEquals(2, buckets[1][2]);
        assertEquals(0, buckets[1][3]);
        assertEquals(0, buckets[1][4]);
        assertEquals(3, buckets[2][0]);
        assertEquals(3, buckets[2][1]);
        assertEquals(1, buckets[2][2]);
        assertEquals(0, buckets[2][3]);
        assertEquals(0, buckets[2][4]);
      }
    } finally {
      assertTrue(FileHelper.deleteAll(dir));
    }
  }

  public void testExec() throws IOException {
    final File dir = FileUtils.createTempDir("cnvproducttask6", "test");
    try {
      final File baseSamFile = new File(dir, "base.sam");
      FileUtils.stringToFile(FileHelper.resourceToString("com/rtg/variant/cnv/resources/babySam.sam"), baseSamFile);
      final File targetSamFile = new File(dir, "target.sam");
      FileUtils.stringToFile(FileHelper.resourceToString("com/rtg/variant/cnv/resources/babySamTarget.sam"), targetSamFile);
      final File output = new File(dir, "output");
      final CnvProductParams params = CnvProductParams.builder().bucketSize(5).mappedBase(Arrays.asList(baseSamFile))
              .mappedTarget(Arrays.asList(targetSamFile)).outputParams(new OutputParams(output, false)).create();
      final CnvProductTask task = new CnvProductTask(params, TestUtils.getNullOutputStream());
      task.run();
      final String bedFile = FileUtils.fileToString(new File(output, "cnv.bed"));
      final String ratioFile = FileUtils.fileToString(new File(output, "cnv.ratio"));
      final String bedBody = stripComments(bedFile);
      final String ratioBody = stripComments(ratioFile);
      final String bedComments = getComments(bedFile);
      final String ratioComments = getComments(ratioFile);
      assertEquals(StringUtils.convertLineEndings(FileHelper.resourceToString("com/rtg/variant/cnv/resources/cnv1.cnv")), bedBody);
      assertEquals(StringUtils.convertLineEndings(FileHelper.resourceToString("com/rtg/variant/cnv/resources/cnv1.ratio")), ratioBody);
      TestUtils.containsAll(bedComments, "#CL\t", "#Version ", ", cnv v",
          "#Seq\tstart\tend\tcnv\tmean",
          "#Seq\tstart\tend\tnregion",
          "#Seq\tstart\tend\tgdelregion",
          "#Seq\tstart\tend\tmodel[0-5] mean");
      TestUtils.containsAll(ratioComments, "#CL\t", "#Version ", ", cnv v",
          "#Seq\tposition\tratio");
    } finally {
      assertTrue(FileHelper.deleteAll(dir));
    }
  }

  public void testDictionaryError() throws IOException {
    final MemoryPrintStream ps = new MemoryPrintStream();
    Diagnostic.setLogStream(ps.printStream());
    try (final TestDirectory dir = new TestDirectory()) {
      final File baseSamFile = new File(dir, "base.sam");
      FileUtils.stringToFile(FileHelper.resourceToString("com/rtg/variant/cnv/resources/babySam.sam"), baseSamFile);
      final File targetSamFile = new File(dir, "target.sam");
      FileUtils.stringToFile(FileHelper.resourceToString("com/rtg/variant/cnv/resources/otherDict.sam"), targetSamFile);
      final File output = new File(dir, "output");
      final CnvProductParams params = CnvProductParams.builder().bucketSize(5).mappedBase(Arrays.asList(baseSamFile))
              .mappedTarget(Arrays.asList(targetSamFile)).outputParams(new OutputParams(output, false)).create();
      final CnvProductTask task = new CnvProductTask(params, TestUtils.getNullOutputStream());
      try {
        task.exec();
        fail();
      } catch (final Exception e) {
        //expected
        assertTrue(ps.toString(), ps.toString().contains("The SAM file \"" + baseSamFile.getPath() + "\" cannot be merged with the SAM file \"" + targetSamFile.getPath() + "\" because their headers are incompatible."));
      }
    } finally {
      Diagnostic.setLogStream();
    }
  }

  private String stripComments(String file) {
    final StringBuilder bd = new StringBuilder();
    final String[] lines = file.split(StringUtils.LS);
    for (String line : lines) {
      if (!line.startsWith("#")) {
        bd.append(line);
        bd.append(StringUtils.LS);
      }
    }
    return bd.toString();
  }

  private String getComments(String file) {
    final StringBuilder bd = new StringBuilder();
    final String[] lines = file.split(StringUtils.LS);
    for (String line : lines) {
      if (line.startsWith("#")) {
        bd.append(line);
        bd.append(StringUtils.LS);
      }
    }
    return bd.toString();
  }

  public void testBloodyBigBucket() throws Exception {
    try (TestDirectory dir = new TestDirectory()) {
      final CnvProductParams params = CnvProductParams.builder().bucketSize(Integer.MAX_VALUE).create();
      try (SamFilterIterator it = new SamFilterIterator(ThreadedMultifileIteratorTest.getIterator(dir, "com/rtg/variant/cnv/resources/babySam.sam"), new DefaultSamFilter(params.filterParams()))) {
        final CnvProductTask task = new CnvProductTask(params, TestUtils.getNullOutputStream());
        task.setSequenceLengths(CnvProductTask.makeTemplateLengths(it.header().getSequenceDictionary()));
        final int[][] buckets = task.chunkSamFile(it);
        assertEquals(3, buckets.length);
        assertEquals(1, buckets[0].length);
        assertEquals(1, buckets[1].length);
        assertEquals(1, buckets[2].length);
      }
    }
  }
}
