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
import java.io.StringWriter;

import com.rtg.launcher.AbstractCli;
import com.rtg.launcher.AbstractCliTest;
import com.rtg.util.io.MemoryPrintStream;
import com.rtg.util.io.TestDirectory;
import com.rtg.util.test.FileHelper;

/**
 */
public class NcbiTaxDumpReaderCliTest extends AbstractCliTest {

  @Override
  protected AbstractCli getCli() {
    return new NcbiTaxDumpReaderCli();
  }

  public void testProcess() throws IOException {
    try (TestDirectory dir = new TestDirectory()) {
      for (String f : new String[]{"nodes.dmp", "names.dmp", "merged.dmp", "delnodes.dmp"}) {
        FileHelper.resourceToFile("com/rtg/taxonomy/resources/" + f, new File(dir, f));
      }
      try (final MemoryPrintStream ps = new MemoryPrintStream()) {
        final NcbiTaxDumpReaderCli.NcbiTaxDumpReader ntdr = new NcbiTaxDumpReaderCli.NcbiTaxDumpReader(dir, ps.printStream());
        final Taxonomy tax = ntdr.getTaxonomy();

        assertTrue(tax.isConsistent());
        assertEquals(14, tax.size());

        final StringWriter sw = new StringWriter();
        tax.write(sw);
        mNano.check("tree2.tsv", sw.toString(), false);
      }
    }
  }

  public void testProcessCli() throws IOException {
    try (TestDirectory dir = new TestDirectory()) {
      for (String f : new String[]{"nodes.dmp", "names.dmp", "merged.dmp", "delnodes.dmp"}) {
        FileHelper.resourceToFile("com/rtg/taxonomy/resources/" + f, new File(dir, f));
      }
      try (final MemoryPrintStream ps = new MemoryPrintStream();
          final MemoryPrintStream psErr = new MemoryPrintStream()) {
        assertEquals(0, getCli().mainInit(new String[] {dir.getPath()}, ps.outputStream(), psErr.printStream()));

        mNano.check("tree2.tsv", ps.toString(), false);
      }
    }
  }
}
