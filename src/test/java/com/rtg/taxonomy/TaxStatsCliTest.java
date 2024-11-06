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
      final String err = checkHandleFlagsErr(tmp.getPath());
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
      FileUtils.stringToFile(TAX, new File(sdf, TaxonomyUtils.TAXONOMY_FILE));
      FileUtils.stringToFile(TAX_LOOK, new File(sdf, TaxonomyUtils.TAXONOMY_TO_SEQUENCE_FILE));
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
      FileUtils.stringToFile(TAX, new File(sdf, TaxonomyUtils.TAXONOMY_FILE));
      FileUtils.stringToFile(TAX_LOOK_DIRE, new File(sdf, TaxonomyUtils.TAXONOMY_TO_SEQUENCE_FILE));
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
      FileUtils.stringToFile(TAX, new File(sdf, TaxonomyUtils.TAXONOMY_FILE));
      FileUtils.stringToFile(TAX_LOOK_DIRE2, new File(sdf, TaxonomyUtils.TAXONOMY_TO_SEQUENCE_FILE));
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
      assertTrue(err.toString(), err.toString().contains("does not seem to contain a valid SDF"));
    }
  }


}
