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
package com.rtg.index.hash.ngs;

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.PrintStream;

import com.rtg.index.hash.ngs.general.Skeleton;
import com.rtg.launcher.ISequenceParams;
import com.rtg.launcher.HashingRegion;
import com.rtg.launcher.SequenceParams;
import com.rtg.mode.SequenceMode;
import com.rtg.reader.ReaderTestUtils;
import com.rtg.util.TestUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.io.FileUtils;
import com.rtg.util.test.FileHelper;

import junit.framework.TestCase;

/**
 * Test for corresponding class
 */
public class ProteinNgsHashLoopTest extends TestCase {

  private static class FakeReadCall implements ReadCall {

    @Override
    public void readCall(final int id, final long hash, final int index) {
    }
  }

  public void testReadLoop() throws Exception {
    Diagnostic.setLogStream();
    final File dir = FileUtils.createTempDir("proteinhashloop", "test");
    try {
      final File readDir = new File(dir, "read");
      final File templateDir = new File(dir, "template");
      final String read = "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa";
      final String template = "KWRKNRKSKKNQRNYNHDAADA";
      final int adjLength = read.length() / 3 - 1;
      ReaderTestUtils.getReaderDNA(">a\n" + read, readDir, null).close();
      ReaderTestUtils.getReaderProtein(">b\n" + template, templateDir).close();
      final ProteinNgsHashLoop loop = new ProteinNgsHashLoop(adjLength, adjLength);
      final FakeProteinMask mask = new FakeProteinMask(new Skeleton(adjLength, adjLength, 0, 0, 1), new FakeReadCall(), new ImplementHashFunctionTest.TemplateCallMock());
      final ISequenceParams readParams = SequenceParams.builder().directory(readDir).mode(SequenceMode.TRANSLATED).create();
      try {
        loop.readLoop(readParams, mask, ReadEncoder.SINGLE_END, false);
      } finally {
        readParams.close();
      }
      assertEquals(6, mask.mReadCalls);
      final ISequenceParams templateParams = SequenceParams.builder().directory(templateDir).mode(SequenceMode.PROTEIN).create();
      try {
        loop.templateLoop(templateParams, mask);
      } finally {
        templateParams.close();
      }
      assertEquals(13, mask.mTemplateCalls);
    } finally {
      assertTrue(FileHelper.deleteAll(dir));
    }
  }

  public void testTemplateLoopMultiCore() throws Exception {
    Diagnostic.setLogStream();
    final File dir = FileUtils.createTempDir("proteinhashloop", "test");
    try {
      final File readDir = new File(dir, "read");
      final File templateDir = new File(dir, "template");
      final String read = "aaatggcgcaaaaacagaaagtcgaaaaaaaatcaa";
      final String template = TestUtils.dnaToProtein("AAATGGCGCAAAAACAGAAAGTCGAAAAAAAATCAAAGAAATTATAACCACGACGCAGCAGACGCAG");
      final int adjLength = read.length() / 3 - 1;
      ReaderTestUtils.getReaderDNA(">b\n" + read, readDir, null).close();
      ReaderTestUtils.getReaderProtein(">b\n" + template + "\n>c\n" + template + "\n>d\n" + template + "\n>e\n" + template, templateDir).close();
      final ProteinNgsHashLoop loop = new ProteinNgsHashLoop(adjLength, adjLength);
      final FakeProteinMask mask = new FakeProteinMask(new Skeleton(adjLength, adjLength, 0, 0, 1), new FakeReadCall(), new ImplementHashFunctionTest.TemplateCallMock());
      final ISequenceParams readParams = SequenceParams.builder().directory(readDir).mode(SequenceMode.TRANSLATED).create();
      try {
        loop.readLoop(readParams, mask, ReadEncoder.SINGLE_END, false);
      } finally {
        readParams.close();
      }
      final ISequenceParams templateParams = SequenceParams.builder().directory(templateDir).mode(SequenceMode.PROTEIN).create();
      try {
        final ByteArrayOutputStream baos = new ByteArrayOutputStream();
        final PrintStream ps = new PrintStream(baos);

        //
        Diagnostic.setLogStream(ps);
        try {
          loop.templateLoopMultiCore(templateParams, mask, 4, 10);
          assertEquals(4, mask.mClones);
          ps.flush();
          final String diagstring = baos.toString();
          assertTrue(diagstring.contains("parent Terminating"));
          assertTrue(diagstring.contains("parent Finished"));
          assertTrue(diagstring.contains("Scheduling"));
          assertTrue(diagstring.contains("Thread Search 1 Start"));
          assertTrue(diagstring.contains(" s"));
        } finally {
          Diagnostic.setLogStream();
        }
      } finally {
        templateParams.close();
      }

      final FakeProteinMask mask2 = new FakeProteinMask(new Skeleton(adjLength, adjLength, 0, 0, 1), new FakeReadCall(), new ImplementHashFunctionTest.TemplateCallMock());
      final ISequenceParams templateParams2 = SequenceParams.builder().directory(templateDir).mode(SequenceMode.PROTEIN).region(new HashingRegion(0, 0)).create();
      try {
        loop.templateLoopMultiCore(templateParams2, mask2, 4, 10);
      } finally {
        templateParams.close();
      }
      assertEquals(0, mask2.mClones);
    } finally {
      assertTrue(FileHelper.deleteAll(dir));
    }
  }
}
