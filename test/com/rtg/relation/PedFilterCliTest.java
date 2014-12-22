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

import com.rtg.launcher.AbstractCli;
import com.rtg.launcher.AbstractCliTest;
import com.rtg.util.StringUtils;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.TestDirectory;

/**
 */
public class PedFilterCliTest extends AbstractCliTest {

  @Override
  protected AbstractCli getCli() {
    return new PedFilterCli();
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

  public void testFile() throws IOException {
    try (final TestDirectory dir = new TestDirectory("pedfilter")) {
      final File relationFile = new File(dir, "relationshipfile.ped");
      final GenomeRelationships ped = makeTestPed();
      FileUtils.stringToFile(PedFileParser.toString(ped), relationFile);

      String output = checkMainInitOk("--vcf", relationFile.toString());
      mNano.check("pedfilter-tovcf.txt", StringUtils.grepMinusV(output, "^##(RUN-ID|CL|fileDate|source)"));
    }
  }

}
