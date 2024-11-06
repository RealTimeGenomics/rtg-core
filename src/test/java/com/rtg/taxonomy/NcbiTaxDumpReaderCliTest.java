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
