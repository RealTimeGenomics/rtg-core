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
