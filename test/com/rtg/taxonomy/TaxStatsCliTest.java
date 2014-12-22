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

package com.rtg.taxonomy;

import java.io.File;
import java.io.IOException;

import com.rtg.launcher.AbstractCli;
import com.rtg.launcher.AbstractCliTest;
import com.rtg.reader.ReaderTestUtils;
import com.rtg.util.TestUtils;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.MemoryPrintStream;
import com.rtg.util.io.TestDirectory;

/**
 */
public class TaxStatsCliTest extends AbstractCliTest {

  @Override
  protected AbstractCli getCli() {
    return new TaxStatsCli();
  }

  public void testValidator() throws IOException {
    try (final TestDirectory dir = new TestDirectory()) {
      final File tmp = new File(dir, "tmp");
      assertTrue(tmp.createNewFile());
      String err = checkHandleFlagsErr(tmp.getPath());
      TestUtils.containsAll(err, "The specified file,", ", is not an SDF.");
      checkHandleFlags(tmp.getPath());
    }
  }

  private static final String REF = ">seq1\nacgt\n>seq2\nacgt\n>seq3\nacgt\n>seq4\nacgt";

  private static final String TAX = "#RTG taxonomy version 1.0\n"
      + "#taxID\tparentID\trank\tname\n"
      + "1\t-1\tno rank\troot\n"
      + "2\t1\tno rank\tbact\n"
      + "3\t2\tno rank\tstrainA\n"
      + "4\t2\tstrain\tstrainB\n"
      + "5\t2\tstrain\tstrainC\n";

  private static final String TAX_LOOK = "2\tseq1\n"
      + "3\tseq2\n"
      + "4\tseq3\n"
      + "5\tseq4\n";

  public void testBasic() throws IOException {
    try (final TestDirectory dir = new TestDirectory()) {
      final File sdf = ReaderTestUtils.getDNADir(REF, new File(dir, "sdf"));
      FileUtils.stringToFile(TAX, new File(sdf, Taxonomy.TAXONOMY_FILE));
      FileUtils.stringToFile(TAX_LOOK, new File(sdf, SequenceToTaxonIds.TAXONOMY_TO_SEQUENCE_FILE));
      final MemoryPrintStream out = new MemoryPrintStream();
      final MemoryPrintStream err = new MemoryPrintStream();
      assertEquals(0, getCli().mainInit(new String[]{sdf.getPath()}, out.outputStream(), err.printStream()));
      mNano.check("taxVerifyErr", err.toString(), true);
      mNano.check("taxVerifyOut", out.toString(), true);
    }
  }

  private static final String TAX_LOOK_DIRE = "2\tseq1\n"
      + "3\tseq2\n"
      + "4\tseq3\n"
      + "5\tseq6\n";

  public void testDire() throws IOException {
    try (final TestDirectory dir = new TestDirectory()) {
      final File sdf = ReaderTestUtils.getDNADir(REF, new File(dir, "sdf"));
      FileUtils.stringToFile(TAX, new File(sdf, Taxonomy.TAXONOMY_FILE));
      FileUtils.stringToFile(TAX_LOOK_DIRE, new File(sdf, SequenceToTaxonIds.TAXONOMY_TO_SEQUENCE_FILE));
      final MemoryPrintStream out = new MemoryPrintStream();
      final MemoryPrintStream err = new MemoryPrintStream();
      assertEquals(1, getCli().mainInit(new String[]{sdf.getPath()}, out.outputStream(), err.printStream()));
      mNano.check("taxVerifyDireErr", err.toString(), true);
    }
  }

  private static final String TAX_LOOK_DIRE2 = "2\tseq1\n"
    + "3\tseq2\n"
    + "4\tseq3\n"
    + "6\tseq4\n";

  public void testDire2() throws IOException {
    try (final TestDirectory dir = new TestDirectory()) {
      final File sdf = ReaderTestUtils.getDNADir(REF, new File(dir, "sdf"));
      FileUtils.stringToFile(TAX, new File(sdf, Taxonomy.TAXONOMY_FILE));
      FileUtils.stringToFile(TAX_LOOK_DIRE2, new File(sdf, SequenceToTaxonIds.TAXONOMY_TO_SEQUENCE_FILE));
      final MemoryPrintStream out = new MemoryPrintStream();
      final MemoryPrintStream err = new MemoryPrintStream();
      assertEquals(1, getCli().mainInit(new String[]{sdf.getPath()}, out.outputStream(), err.printStream()));
      mNano.check("taxVerifyDire2Err", err.toString(), true);
    }
  }

  public void testNotSdf() throws IOException {
    try (final TestDirectory dir = new TestDirectory()) {
      final MemoryPrintStream out = new MemoryPrintStream();
      final MemoryPrintStream err = new MemoryPrintStream();
      assertEquals(1, getCli().mainInit(new String[]{dir.getPath()}, out.outputStream(), err.printStream()));
      assertTrue(err.toString().contains("is not an SDF"));
    }
  }


}
