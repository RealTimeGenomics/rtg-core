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

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.Arrays;

import com.rtg.launcher.SequenceParams;
import com.rtg.ngs.blocking.MapQScoringReadBlocker;
import com.rtg.ngs.blocking.MapQScoringReadBlockerSynch;
import com.rtg.reader.Names;
import com.rtg.reader.PrereadArm;
import com.rtg.reader.ReadHelper;
import com.rtg.reader.ReaderTestUtils;
import com.rtg.reader.SequencesReader;
import com.rtg.reference.Sex;
import com.rtg.util.intervals.LongRange;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.NoTalkbackSlimException;
import com.rtg.util.io.FileUtils;
import com.rtg.util.test.FileHelper;

import junit.framework.TestCase;

/**
 * Test class
 */
public class SharedResourcesTest extends TestCase {

  /**
   */
  public SharedResourcesTest(final String name) {
    super(name);
  }

  public final void testSharedResources() throws IOException {
    Diagnostic.setLogStream();
    final File topLevel = FileUtils.createTempDir("sharedresource", "test");
    try {
      final File leftFile = new File(topLevel, "left");
      final SequencesReader left = ReaderTestUtils.getReaderDNA(">testLeft\nacgt", leftFile, null, true);
      final File rightFile = new File(topLevel, "right");
      final SequencesReader right = ReaderTestUtils.getReaderDNA(">testRight\nacct", rightFile, null, true);
      final File refFile = new File(topLevel, "ref");
      final SequencesReader ref = ReaderTestUtils.getReaderDNA(">testRef\nacct", refFile, null, true);

      final MapQScoringReadBlocker srb = new MapQScoringReadBlocker(1, 2);
      final Names names = new Names(refFile, LongRange.NONE);
      final SharedResources sr = new SharedResources(left, right, ref, srb, names, SharedResources.createHeader(ref, names, left, null, Sex.EITHER, true), null);
      assertTrue(srb == sr.getBlocker());
      assertTrue(sr.isPairedEnd());
      try (SequencesReader rightReaderCopy = sr.secondReaderCopy(); SequencesReader leftReaderCopy = sr.firstReaderCopy(); SequencesReader refReaderCopy = sr.templateReaderCopy()) {
        assertTrue(Arrays.equals(new byte[]{1, 2, 3, 4}, ReadHelper.getRead(leftReaderCopy, 0)));

        assertTrue(Arrays.equals(new byte[]{1, 2, 2, 4}, ReadHelper.getRead(rightReaderCopy, 0)));
        assertNull(ReadHelper.getQual(leftReaderCopy, 0));
        assertNull(ReadHelper.getQual(rightReaderCopy, 0));
        assertNull(ReadHelper.getQual(refReaderCopy, 0));
        assertTrue(names == sr.names());
      }
    } finally {
      assertTrue(FileHelper.deleteAll(topLevel));
    }
  }

  public final void testSharedResourcesClose() throws IOException {
    final SharedResources sr = new SharedResources(null, null, null, null, null, null, null);
    sr.close();
  }

  public void testNulls() {
    final SharedResources r = new SharedResources(null, null, null, null, null, null, null);
    assertNull(r.getBlocker());
    assertNull(r.firstReaderCopy());
    assertNull(r.secondReaderCopy());
    assertFalse(r.isPairedEnd());
  }

  public void testGeneratingFromParams() throws IOException {
    Diagnostic.setLogStream();
    final File topLevel = FileUtils.createTempDir("sharedresource", "test");
    try {
      final File leftFile = new File(topLevel, "left");
      final StringBuilder elevenReads = new StringBuilder();
      for (int i = 0; i < 11; ++i) {
        elevenReads.append(">testLeft-").append(i).append("\nacgt\n");
      }
      ReaderTestUtils.getReaderDNA(elevenReads.toString(), leftFile, null);
      final File rightFile = new File(topLevel, "right");
      ReaderTestUtils.getReaderDNA(">testRight\nacct", rightFile, null);

      final NgsFilterParams filterParams = new NgsFilterParams.NgsFilterParamsBuilder().maxTopResults(10).create();
      final NgsOutputParams outputParams = new NgsOutputParamsBuilder().filterParams(filterParams).create();
      final NgsParams p = new NgsParamsBuilder().buildFirstParams(new SequenceParams.SequenceParamsBuilder().directory(leftFile).create())
                                                          .buildSecondParams(new SequenceParams.SequenceParamsBuilder().directory(rightFile).create())
                                                          .searchParams(new SequenceParams.SequenceParamsBuilder().directory(leftFile).create())
                                                          .numberThreads(2)
                                                          .outputParams(outputParams)
                                                          .create();

      final SharedResources sr = SharedResources.generateSharedResources(p);
      try (SequencesReader left = sr.firstReaderCopy(); SequencesReader right = sr.secondReaderCopy()) {
        assertTrue(Arrays.equals(new byte[]{1, 2, 3, 4}, ReadHelper.getRead(left, 0)));
        assertTrue(Arrays.equals(new byte[]{1, 2, 2, 4}, ReadHelper.getRead(right, 0)));
        assertNull(ReadHelper.getQual(left, 0));
        assertNull(ReadHelper.getQual(right, 0));
        assertNotNull(sr.getBlocker());
        assertTrue(sr.getBlocker() instanceof MapQScoringReadBlockerSynch);
        assertNotNull(sr.names());

        final NgsParams p2 = new NgsParamsBuilder().buildFirstParams(new SequenceParams.SequenceParamsBuilder().directory(leftFile).create())
          .buildSecondParams(new SequenceParams.SequenceParamsBuilder().directory(rightFile).create())
          .searchParams(new SequenceParams.SequenceParamsBuilder().directory(leftFile).create())
          .numberThreads(1)
          .outputParams(outputParams)
          .create();
        final SharedResources sr2 = SharedResources.generateSharedResources(p2);
        try {
          assertTrue(sr2.getBlocker() instanceof MapQScoringReadBlockerSynch);
        } finally {
          sr2.close();
        }
      } finally {
        sr.close();


      }
    } finally {

      assertTrue(FileHelper.deleteAll(topLevel));
    }
  }

  public void testValidation() throws Exception {
    final File tmp = FileUtils.createTempDir("sharedresources", "blkjer");
    final ByteArrayOutputStream os = new ByteArrayOutputStream();
    final PrintStream ps = new PrintStream(os);
    Diagnostic.setLogStream(ps);
    try {
      try (SequencesReader reader = ReaderTestUtils.getReaderDNAFastqCG("@test0\nacgtacgtacgtacgtacgtacgtacgtacgtacg\n+\n###################################\n@test1\na\n+\n#\n", tmp, PrereadArm.LEFT)) {
        new SharedResources(reader, reader, reader, null, null, null, null);
        fail();
      } catch (final NoTalkbackSlimException ntse) {
        ntse.logException();
        ps.flush();
        assertTrue(os.toString().contains("Sequences of varying length found in: \"" + tmp.toString() + "\""));
      }

    } finally {
      Diagnostic.setLogStream();
      FileHelper.deleteAll(tmp);
    }
  }

}
