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
package com.rtg.variant.coverage;

import static com.rtg.util.StringUtils.LS;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.util.ArrayList;
import java.util.List;

import com.rtg.launcher.OutputParams;
import com.rtg.launcher.ReaderParams;
import com.rtg.launcher.SequenceParams;
import com.rtg.mode.SequenceMode;
import com.rtg.reader.ReaderTestUtils;
import com.rtg.sam.SamFilterParams;
import com.rtg.sam.SharedSamConstants;
import com.rtg.util.TestUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.io.FileUtils;
import com.rtg.util.test.FileHelper;
import com.rtg.util.test.params.TestParams;
import com.rtg.variant.coverage.CoverageParams.CoverageParamsBuilder;

import junit.framework.TestCase;

/**
 */
public class CoverageParamsTest extends TestCase {

  protected File mDir;

  @Override
  public void setUp() throws IOException {
    Diagnostic.setLogStream();
    mDir = FileHelper.createTempDirectory();
  }

  @Override
  public void tearDown() {
    assertTrue(FileHelper.deleteAll(mDir));
    mDir = null;
  }

  private static final String TEST_OUTPUT = "coveragetestoutput";

  private ReaderParams makeGenome() throws IOException {
    final File subjectsDir = FileUtils.createTempDir("test", "coverageparams", mDir);
    ReaderTestUtils.getReaderDNA(">t\nacgt", subjectsDir, null).close();
    return SequenceParams.builder().directory(subjectsDir).mode(SequenceMode.UNIDIRECTIONAL).create().readerParams();
  }

  public void testOmnes() {
    new TestParams(CoverageParams.class, CoverageParamsBuilder.class).check();
  }

  public void testDefaultParams() throws IOException {
    final File outDir = new File(mDir, "output");
    final CoverageParams cp = CoverageParams.builder().outputParams(new OutputParams(outDir, false)).genome(makeGenome()).create();
    assertFalse(cp.tsvOutput());
    assertTrue(cp.bedOutput());
    assertFalse(cp.blockCompressed());
    assertFalse(cp.errorRates());
    assertEquals(0, cp.smoothing());
    assertEquals(1, cp.minimumCoverageThreshold());
    assertEquals(outDir, cp.directory());
    assertEquals(cp.outFile(), cp.file("coverage.bed"));
    assertEquals("coverage.bed", cp.outFile().getName());
    FileUtils.ensureOutputDirectory(cp.directory());
    final OutputStream out = cp.bedStream();
    out.write("test".getBytes());
    out.close();
    assertTrue(cp.outFile().exists());
    assertEquals("test", FileUtils.fileToString(cp.outFile()));
  }

  public void testTsv() throws IOException {
    final File outDir = new File(mDir, "output");
    final CoverageParams cp = CoverageParams.builder().outputParams(new OutputParams(outDir, false)).genome(makeGenome()).tsvOutput(true).create();
    assertTrue(cp.tsvOutput());
    assertEquals("coverage.tsv", cp.outFile().getName());
    assertFalse(cp.bedOutput());
  }

  private static final String[] BASE_PARAMS = {
    "CoverageParams mapped reads=",
    "smoothing=0",
    "error rates=" + Boolean.FALSE.toString(),
    "    ReaderParams directory=",
    "OutputParams output directory=",
    "zip=" + Boolean.TRUE + LS,
    "min coverage threshold=3"
  };

  CoverageParams getCoverageParams(final File outFile, final List<File> mapped) throws IOException {
    return CoverageParams.builder().outputParams(new OutputParams(outFile, true))
                                   .mapped(mapped).genome(makeGenome()).smoothing(0)
                                   .ioThreads(1).minimumCoverageThreshold(3).tsvOutput(false)
                                   .name("CoverageParams").errorRates(false)
                                   .filterParams(SamFilterParams.builder().create())
                                   .create();
  }

  public void testOk1() throws Exception {
    final File map = File.createTempFile("testok1", "coverageParams");
    try {
      FileUtils.stringToFile(SharedSamConstants.SAM9, map);
      assertTrue(map.isFile());

      final CoverageParams ccp;
      final File outFile = new File(mDir, TEST_OUTPUT);
      assertTrue(outFile.mkdir());

      final List<File> mapped = new ArrayList<>();
      mapped.add(map);

      ccp = getCoverageParams(outFile, mapped);

      assertFalse(ccp.errorRates());
      assertNotNull(ccp.filterParams());
      assertEquals(1, ccp.ioThreads());
      assertEquals(3, ccp.minimumCoverageThreshold());
      final String ccs = ccp.toString();
      //System.err.println(ccs);
      TestUtils.containsAll(ccs, BASE_PARAMS);
      ccp.close();
    } finally {
      assertTrue(map.delete());
    }
  }

}
