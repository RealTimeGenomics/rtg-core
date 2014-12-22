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
package com.rtg.reader;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

import com.rtg.util.io.FileUtils;
import com.rtg.util.test.FileHelper;

import junit.framework.TestCase;

/**
 */
public class SourceTemplateReadWriterTest extends TestCase {

  public void test() throws IOException {
    final File temp = FileUtils.createTempDir("sourcetemplate", "tempdir");
    try {
      final File templateFile = new File(temp, "simTemplates");
      SourceTemplateReadWriter.writeTemplateMappingFile(temp, null);
      assertFalse(templateFile.exists());
      SdfId[] resultArray = SourceTemplateReadWriter.readTemplateMap(temp);
      assertNull(resultArray);
      final SdfId[] testArray = {new SdfId(Long.MIN_VALUE), new SdfId(0), new SdfId(Long.MAX_VALUE)};
      SourceTemplateReadWriter.writeTemplateMappingFile(temp, testArray);
      assertTrue(templateFile.exists());
      resultArray = SourceTemplateReadWriter.readTemplateMap(temp);
      assertEquals(testArray.length, resultArray.length);
      for (int i = 0; i < testArray.length; i++) {
        assertEquals(testArray[i], resultArray[i]);
      }
      assertTrue(templateFile.delete());
      assertTrue(templateFile.createNewFile());
      try (BufferedWriter bw = new BufferedWriter(new FileWriter(templateFile))) {
        bw.write("blah");
      }
      try {
        SourceTemplateReadWriter.readTemplateMap(temp);
        fail();
      } catch (CorruptSdfException e) {
        assertEquals("Malformed simulator template map line blah", e.getMessage());
      }
    } finally {
      FileHelper.deleteAll(temp);
    }
  }

  public void test2() throws IOException {
    final File temp = FileUtils.createTempDir("sourcetemplate", "tempdir");
    try {
      final File templateFile = new File(temp, "mutTemplate");
      SourceTemplateReadWriter.writeMutationMappingFile(temp, null);
      assertFalse(templateFile.exists());
      SdfId result = SourceTemplateReadWriter.readMutationMap(temp);
      assertNull(result);
      final SdfId testSdf = new SdfId(Long.MIN_VALUE);
      SourceTemplateReadWriter.writeMutationMappingFile(temp, testSdf);
      assertTrue(templateFile.exists());
      result = SourceTemplateReadWriter.readMutationMap(temp);
      assertEquals(testSdf, result);
      final File newDir = new File(temp, "newDir");
      assertTrue(newDir.mkdir());
      SourceTemplateReadWriter.copyMutationMappingFile(temp, newDir);
      result = SourceTemplateReadWriter.readMutationMap(newDir);
      assertEquals(testSdf, result);
      assertTrue(templateFile.delete());
      assertTrue(templateFile.createNewFile());
      try (BufferedWriter bw = new BufferedWriter(new FileWriter(templateFile))) {
        bw.write("blah");
      }
      try {
        SourceTemplateReadWriter.readMutationMap(temp);
        fail();
      } catch (CorruptSdfException e) {
        assertEquals("Malformed simulator template map line blah", e.getMessage());
      }
    } finally {
      FileHelper.deleteAll(temp);
    }
  }
}
