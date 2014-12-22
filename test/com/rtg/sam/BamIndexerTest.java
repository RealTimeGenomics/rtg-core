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
package com.rtg.sam;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.InputStream;

import com.rtg.tabix.IndexTestUtils;
import com.rtg.util.Resources;

import junit.framework.TestCase;

/**
 *
 */
public class BamIndexerTest extends TestCase {

  public void test() throws Exception {
    final ByteArrayOutputStream os = new ByteArrayOutputStream();
    try {
      try (InputStream is = Resources.getResourceAsStream("com/rtg/sam/resources/bam.bam")) {
        BamIndexer.saveBamIndex(is, os);
      }
    } finally {
      os.close();
    }
    final String myBai = IndexTestUtils.bamIndexToUniqueString(new ByteArrayInputStream(os.toByteArray()));
    final String exp;
    try (InputStream baiIs = Resources.getResourceAsStream("com/rtg/sam/resources/bam.bam.bai")) {
      exp = IndexTestUtils.bamIndexToUniqueString(baiIs);
    }
    assertEquals(exp, myBai);
  }

  public void test2() throws Exception {
    final ByteArrayOutputStream os = new ByteArrayOutputStream();
    try {
      try (InputStream is = Resources.getResourceAsStream("com/rtg/sam/resources/mated.bam")) {
        BamIndexer.saveBamIndex(is, os);
      }
    } finally {
      os.close();
    }
    final String myBai = IndexTestUtils.bamIndexToUniqueString(new ByteArrayInputStream(os.toByteArray()));
    final String exp;
    try (InputStream baiIs = Resources.getResourceAsStream("com/rtg/sam/resources/mated.bam.bai")) {
      exp = IndexTestUtils.bamIndexToUniqueString(baiIs);
    }
    assertEquals(exp, myBai);
  }

  public void test3() throws Exception {
    final ByteArrayOutputStream os = new ByteArrayOutputStream();
    try {
      try (InputStream is = Resources.getResourceAsStream("com/rtg/sam/resources/multiSequence.bam")) {
        BamIndexer.saveBamIndex(is, os);
      }
    } finally {
      os.close();
    }
    final String myBai = IndexTestUtils.bamIndexToUniqueString(new ByteArrayInputStream(os.toByteArray()));
    final String exp;
    try (InputStream baiIs = Resources.getResourceAsStream("com/rtg/sam/resources/multiSequence.bam.bai")) {
      exp = IndexTestUtils.bamIndexToUniqueString(baiIs);
    }
    assertEquals(exp, myBai);
  }

  //test combined mapped and unmapped file
  public void testMixed() throws Exception {
    final ByteArrayOutputStream os = new ByteArrayOutputStream();
    try {
      try (InputStream is = Resources.getResourceAsStream("com/rtg/sam/resources/mixed.bam")) {
        BamIndexer.saveBamIndex(is, os);
      }
    } finally {
      os.close();
    }
    final String myBai = IndexTestUtils.bamIndexToUniqueString(new ByteArrayInputStream(os.toByteArray()));
    final String exp;
    try (InputStream baiIs = Resources.getResourceAsStream("com/rtg/sam/resources/mixed.bam.bai")) {
      exp = IndexTestUtils.bamIndexToUniqueString(baiIs);
    }
    assertEquals(exp, myBai);
  }

  public void testIndexFileName() {
    final File f1 = new File("test.bam");
    assertEquals("test.bam.bai", BamIndexer.indexFileName(f1).getName());
    assertEquals("test.bai", BamIndexer.secondaryIndexFileName(f1).getName());
    final File f2 = new File("test");
    assertEquals("test.bai", BamIndexer.indexFileName(f2).getName());
    assertEquals("test.bai", BamIndexer.secondaryIndexFileName(f2).getName());
  }
}
