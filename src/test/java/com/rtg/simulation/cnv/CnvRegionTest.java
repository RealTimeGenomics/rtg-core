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
package com.rtg.simulation.cnv;

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.io.StringReader;

import com.rtg.reader.PrereadType;
import com.rtg.reader.ReaderTestUtils;
import com.rtg.reader.SdfWriter;
import com.rtg.reader.SequencesReader;
import com.rtg.reader.SequencesReaderFactory;
import com.rtg.simulation.CnvFileChecker;
import com.rtg.util.Constants;
import com.rtg.util.InvalidParamsException;
import com.rtg.util.StringUtils;
import com.rtg.util.TestUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.io.FileUtils;
import com.rtg.util.test.FileHelper;
import com.rtg.util.test.NotRandomRandom;

import junit.framework.TestCase;


/**
 */
public class CnvRegionTest extends TestCase {

  private File mDir;
  private CnvPriorParams mPriors;
  //private NotRandomRandom mRandom;
  private static final String TB = "\t";


  /**
   */
  public CnvRegionTest(String name) {
    super(name);
  }

  public void test() {
    final CnvRegion region = new CnvRegion(0, 40, 100);

    final CnvRegion region2 = new CnvRegion(0, 40, 0);
    region2.mNumCopies = 9;
    assertEquals(0, region2.getCN());

    region.addCopy(region2);
    //System.err.println(region.toString());
    assertEquals("seq: 0 cnved: false start: 40 length: 100 del-1: false del-2: false copies: 0 "
        + "num copies added: 1 num flags added: 1" + StringUtils.LS
        + " heterolist: false",
        region.toString());

    final NotRandomRandom mRandom = new NotRandomRandom();
    final CnvRegion regionPriors = new CnvRegion(0, 40, 100, mPriors);
    assertFalse(regionPriors.isUnpickable(new int[]{1000}, true, -1));
    regionPriors.initializeAsCnved(mRandom);
    assertTrue(regionPriors.isUnpickable(new int[]{1000}, true, -1));
    TestUtils.containsAll(regionPriors.toString(), "copies: 0", "del-1: true", "del-2: true");
    final CnvRegion regionPriors2 = new CnvRegion(0, 40, 100, mPriors);
    regionPriors2.initializeAsCnved(mRandom);
    TestUtils.containsAll(regionPriors2.toString(), "copies: 0", "del-1: true", "del-2: false");
    final CnvRegion regionPriors3 = new CnvRegion(0, 40, 100, mPriors);
    regionPriors3.initializeAsCnved(mRandom);
    TestUtils.containsAll(regionPriors3.toString(), "copies: 1", "del-1: false", "del-2: false");
    final CnvRegion regionPriors6 = new CnvRegion(0, 40, 100, mPriors);
    regionPriors6.initializeAsCnved(mRandom);
    TestUtils.containsAll(regionPriors6.toString(), "copies: 2", "del-1: false", "del-2: false");
    assertEquals(2, CnvRegion.generateNumCopies(mRandom, mPriors, 5));
    assertEquals(4, CnvRegion.generateNumCopies(mRandom, mPriors, 5));
    assertEquals(8, CnvRegion.generateNumCopies(mRandom, mPriors, 5));
    assertEquals(9, CnvRegion.generateNumCopies(mRandom, mPriors, 5));


    //System.err.println(regionPriors.toString());
    final CnvRegion region3 = new CnvRegion(0, 40, 0, 0);
    TestUtils.containsAll(region3.toString(), "del-1: true");
    TestUtils.containsAll(region3.toString(), "del-1: true", "del-2: true");
  }

  public void testGenerate() throws IOException {
    Diagnostic.setLogStream();
    //final File temp = FileHelper.createTempDirectory();
    final String genome = ""
      + ">seq1\n"
      + "AAAAATTTTTGGGGGAAAAATTTTTGGGGGAAAAATTTTTGGGGGAAAAATTTTTGGGGGAAAAA" + StringUtils.LS;
    final File in = ReaderTestUtils.getDNASubDir(genome, mDir);
    try {
      try (SequencesReader dsr = SequencesReaderFactory.createDefaultSequencesReader(in)) {
        final File cnvs = new File(mDir, "test.cnv");
        final File outputDirectory = new File(mDir, "out");
        final File twinDirectory = new File(mDir, "twin");

        try (SdfWriter output = new SdfWriter(outputDirectory, Constants.MAX_FILE_SIZE,
          PrereadType.UNKNOWN, false, true, false, dsr.type())) {
          try (SdfWriter twin = new SdfWriter(twinDirectory, Constants.MAX_FILE_SIZE,
            PrereadType.UNKNOWN, false, true, false, dsr.type())) {
            final CnvSimulator cs = new CnvSimulator(dsr, output, twin, new FileOutputStream(cnvs), new NotRandomRandom(), mPriors, 10, Integer.MAX_VALUE);
            cs.generate();
            final String csstr = cs.toString();
            TestUtils.containsAll(csstr, "CnvSimulator", "No. breakpoints ",
              "Sequence 0 No. regions 2", "priors set");
            final File psFile = new File(mDir, "console");
            final FileOutputStream ps = new FileOutputStream(psFile);
            try {
              cs.outputHistograms(ps);
            } finally {
              ps.flush();
              ps.close();
            }
          }
        }
        final String cnvfilestr = FileUtils.fileToString(cnvs);
        //System.err.println(cnvfilestr);
        TestUtils.containsAll(cnvfilestr, "#Seq" + TB + "start" + TB + "end" + TB + "label" + TB + "cn" + TB + "bp-cn" + TB + "error",
          "seq1" + TB + "0" + TB + "65" + TB + "cnv" + TB + "2" + TB + "0" + TB + "0.0");

        final ByteArrayOutputStream bos = new ByteArrayOutputStream();

        new CnvFileChecker(new PrintStream(bos), new StringReader(cnvfilestr)).check();

        bos.close();

        assertEquals("", bos.toString());

        final int totalLength = CnvSimulatorTest.calculateTotalLength(cnvfilestr);
        try (SequencesReader outputReader = SequencesReaderFactory.createDefaultSequencesReader(outputDirectory)) {
          try (SequencesReader twinReader = SequencesReaderFactory.createDefaultSequencesReader(twinDirectory)) {
            assertEquals(totalLength, outputReader.totalLength() + twinReader.totalLength());
          }
        }
      }
    } finally {
    }
  }
  @Override
  public void setUp() throws IOException, InvalidParamsException {
    mDir = FileHelper.createTempDirectory();
    //mRandom = new NotRandomRandom();
    mPriors = CnvPriorParams.builder().cnvpriors("cnv-default").create();
  }

  @Override
  public void tearDown()  {
    assertTrue(FileHelper.deleteAll(mDir));
    mDir = null;
    //mRandom = null;
    mPriors = null;
  }
}
