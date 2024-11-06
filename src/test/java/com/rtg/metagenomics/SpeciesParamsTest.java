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
package com.rtg.metagenomics;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

import com.rtg.launcher.OutputParams;
import com.rtg.launcher.ReaderParams;
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
      final ReaderParams genomes = SequenceParams.builder().directory(genomeDir).mode(SequenceMode.UNIDIRECTIONAL).create().readerParams();
      final File map = File.createTempFile("testok1", "speciesParams");
      try {
        FileUtils.stringToFile(SharedSamConstants.SAM9, map);
        assertTrue(map.isFile());

        final SpeciesParams ccp;
        try {
          final File outFile = new File(TEST_OUTPUT);
          assertTrue(outFile.mkdir());

          final List<File> mapped = new ArrayList<>();
          mapped.add(map);

          final SpeciesParamsBuilder build = SpeciesParams.builder().mapped(mapped).genome(genomes).referenceMap(null);
          ccp = build.filterParams(SamFilterParams.builder().create()).minIter(10).verbose(false).outputParams(new OutputParams(genomeDir, false)).printAll(true).create();
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
