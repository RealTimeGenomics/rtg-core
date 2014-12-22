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
package com.rtg.metagenomics;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

import com.rtg.launcher.OutputParams;
import com.rtg.launcher.SequenceParams;
import com.rtg.metagenomics.SpeciesParams.SpeciesParamsBuilder;
import com.rtg.mode.SequenceMode;
import com.rtg.reader.ReaderTestUtils;
import com.rtg.sam.SamFilterParams;
import com.rtg.sam.SharedSamConstants;
import com.rtg.util.TestUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.io.FileUtils;
import com.rtg.util.test.FileHelper;
import com.rtg.util.test.params.TestParams;

import junit.framework.TestCase;

/**
 */
public class SpeciesParamsTest extends TestCase {

  private static final String TEST_OUTPUT = "speciestestoutput";

  private static final String[] BASE_PARAMS = {
    "SpeciesParams",
    "SamFilterParams",
  };

  public void testOmnes() {
    new TestParams(SpeciesParams.class, SpeciesParamsBuilder.class).check();
  }

  public void testOk1() throws Exception {
    Diagnostic.setLogStream();
    final File genomeDir = FileUtils.createTempDir("test", "coverageparams");
    try {
      ReaderTestUtils.getReaderDNA(">t\nacgt", genomeDir, null).close();
      final SequenceParams genomes = SequenceParams.builder().directory(genomeDir).mode(SequenceMode.UNIDIRECTIONAL).create();
      final File map = File.createTempFile("testok1", "speciesParams");
      try {
        FileUtils.stringToFile(SharedSamConstants.SAM9, map);
        assertTrue(map.isFile());

        SpeciesParams ccp;
        try {
          final File outFile = new File(TEST_OUTPUT);
          assertTrue(outFile.mkdir());

          final List<File> mapped = new ArrayList<>();
          mapped.add(map);

          final SpeciesParamsBuilder build = SpeciesParams.builder().mapped(mapped).genome(genomes).referenceMap(null);
          ccp = build.filterParams(SamFilterParams.builder().create()).minIter(10).verbose(false).outputParams(new OutputParams(genomeDir, false, false)).printAll(true).create();
          assertNotNull(ccp.filterParams());
          assertEquals(genomes, ccp.genome());
          assertNull(ccp.referenceMap());
          assertFalse(ccp.verbose());
          assertEquals(10, ccp.minIter());
          assertEquals(mapped, ccp.mapped());
          assertNotNull(ccp.filterParams());
          assertEquals(genomeDir, ccp.directory());
          assertEquals(new File(genomeDir, "test"), ccp.file("test"));
          ccp.speciesStream().close();
          assertTrue(new File(genomeDir, "species.tsv").isFile());
          assertEquals(1, ccp.ioThreads());
          assertTrue(ccp.printAll());
          final String ccs = ccp.toString();
          //System.err.println(ccs);
          TestUtils.containsAll(ccs, BASE_PARAMS);
        } finally {
          FileHelper.deleteAll(new File(TEST_OUTPUT));
        }

      } finally {
        assertTrue(map.delete());
      }
    } finally {
      assertTrue(FileHelper.deleteAll(genomeDir));
    }
  }

}
