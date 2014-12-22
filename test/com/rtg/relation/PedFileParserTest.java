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
package com.rtg.relation;

import java.io.File;
import java.io.IOException;

import com.rtg.util.io.FileUtils;
import com.rtg.util.test.FileHelper;
import com.rtg.util.test.NanoRegression;

import junit.framework.TestCase;

/**
 */
public class PedFileParserTest extends TestCase {

  private NanoRegression mNano = null;

  @Override
  public void setUp() {
    mNano = new NanoRegression(this.getClass());
  }
  @Override
  public void tearDown() throws Exception {
    try {
      mNano.finish();
    } finally {
      mNano = null;
    }
  }

  private GenomeRelationships makeTestPed() {
    GenomeRelationships genomeRelationships = new GenomeRelationships();
    genomeRelationships.addGenome("father", GenomeRelationships.SEX_MALE).setProperty(GenomeRelationships.DISEASE_PROPERTY, "true");
    genomeRelationships.addGenome("mother", GenomeRelationships.SEX_FEMALE);
    genomeRelationships.addGenome("child", GenomeRelationships.SEX_MALE).setProperty(GenomeRelationships.DISEASE_PROPERTY, "true");
    genomeRelationships.addParentChild("father", "child");
    genomeRelationships.addParentChild("mother", "child");
    return genomeRelationships;
  }

  public void testToString() throws IOException {
    final GenomeRelationships ped = makeTestPed();
    mNano.check("pednormal", PedFileParser.toString(ped), true);
  }

  public void testFile() throws IOException {
    final File dir = FileUtils.createTempDir("test", "relationshipfile");
    try {
      final File relationFile = new File(dir, "relationshipfile.ped");
      final GenomeRelationships ped = makeTestPed();
      FileUtils.stringToFile(PedFileParser.toString(ped), relationFile);
      final GenomeRelationships ped2 = PedFileParser.loadFile(relationFile);
      mNano.check("pedfile", PedFileParser.toString(ped2));
    } finally {
      assertTrue(FileHelper.deleteAll(dir));
    }
  }

}
