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

package com.rtg.visualization;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.HashMap;

import com.rtg.reader.ReaderTestUtils;
import com.rtg.reader.SdfId;
import com.rtg.util.StringUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.NoTalkbackSlimException;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.MemoryPrintStream;
import com.rtg.util.test.FileHelper;

import junit.framework.TestCase;

/**
 */
public class ReferenceHelperTest extends TestCase {

  /**
   * Test method for {@link com.rtg.visualization.ReferenceHelper#getReference(java.util.HashMap, java.lang.String, int, int)}.
   */
  public final void testGetReference() {
    final HashMap<String, byte[]> map = new HashMap<>();
    map.put("test", new byte[] {1, 2, 3, 4});
    assertEquals("cg", ReferenceHelper.getReference(map, "test", 1, 3));
  }

  /**
   * Test method for {@link com.rtg.visualization.ReferenceHelper#getReferenceFromBytes(byte[], int, int)}.
   */
  public final void testGetReferenceFromBytes() {
    assertEquals("na", ReferenceHelper.getReferenceFromBytes(new byte[] {0, 1, 2, 3}, 0, 2));
  }

  public final void testReferenceFromBytesInclusive() {
    assertEquals("nac", ReferenceHelper.referenceStringFromBytes(new byte[] {0, 1, 2, 3}, 0, 2));
  }


  /**
   * Test method for {@link com.rtg.visualization.ReferenceHelper#loadTemplate(java.io.File, java.util.HashMap, java.io.PrintStream)}.
   */
  public final void testLoadTemplateFile2() throws IOException {
    Diagnostic.setLogStream();
    final File f = FileUtils.createTempDir("loadtemplate", "commonserverutils");
    try {
      final File template = new File(f, "read");
      ReaderTestUtils.getReaderDNA(">name" + StringUtils.LS + "ACGTN", template, new SdfId(1L)).close();
      final HashMap<String, byte[]> m = new HashMap<>();
      final MemoryPrintStream ps = new MemoryPrintStream();
      ReferenceHelper.loadTemplate(template, m, ps.printStream());
      assertEquals("Template loaded : 1" + StringUtils.LS, ps.toString());
      ps.close();
      assertTrue(m.size() == 1);
      final byte[] d = m.get("name");
      assertTrue(Arrays.equals(new byte[] {1, 2, 3, 4, 0}, d));
    } finally {
      assertTrue(FileHelper.deleteAll(f));
    }
  }

  /**
   * Test method for {@link com.rtg.visualization.ReferenceHelper#loadAllTemplate(java.io.File)}.
   */
  public final void testLoadTemplateFile() throws IOException {
    Diagnostic.setLogStream();
    final File f = FileUtils.createTempDir("loadtemplate", "commonserverutils");
    try {
      final File template = new File(f, "read");
      ReaderTestUtils.getReaderDNA(">name" + StringUtils.LS + "ACGTN", template, new SdfId(1L)).close();
      final HashMap<String, byte[]> m = ReferenceHelper.loadAllTemplate(template);
      assertTrue(m.size() == 1);
      final byte[] d = m.get("name");
      assertTrue(Arrays.equals(new byte[] {1, 2, 3, 4, 0}, d));
    } finally {
      assertTrue(FileHelper.deleteAll(f));
    }
  }

  public void testSingleLoad() throws IOException {
    Diagnostic.setLogStream();
    final File f = FileUtils.createTempDir("loadtemplate", "commonserverutils");
    try {
      final File template = new File(f, "read");
      ReaderTestUtils.getReaderDNA(">name" + StringUtils.LS + "ACGTN" + StringUtils.LS
          + ">name2" + StringUtils.LS + "TTTN", template, new SdfId(1L)).close();

      final byte[] d = ReferenceHelper.loadSingleTemplate(template, "name");
      assertTrue(Arrays.equals(new byte[] {1, 2, 3, 4, 0}, d));

      try {
        ReferenceHelper.loadSingleTemplate(template, "wrongname");
      } catch (NoTalkbackSlimException ex) {
        assertTrue(ex.getMessage().contains("Given sequence name not found : wrongname"));
      }
    } finally {
      assertTrue(FileHelper.deleteAll(f));
    }
  }

}
