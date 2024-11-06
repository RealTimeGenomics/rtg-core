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
package com.rtg.variant.bayes.multisample.population;

import java.io.File;

import com.rtg.tabix.TabixIndexer;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.io.MemoryPrintStream;
import com.rtg.util.io.TestDirectory;
import com.rtg.util.test.FileHelper;

import junit.framework.TestCase;

/**
 */
public class AlleleCountsFileReaderTest extends TestCase {

  public void testTextFormat() throws Exception {
    try (TestDirectory td = new TestDirectory()) {
      final String singleLine = "blah\t1\tA\t5\tG\t2\n";

      final File out = new File(td, "snps.ac.gz");
      FileHelper.stringToGzFile(singleLine, out);

      new TabixIndexer(out, new File(td, "snps.ac.gz.tbi")).saveAlleleCountsIndex();

      try (AlleleCountsFileReader afr = AlleleCountsFileReader.openAlleleCountReader(out, null)) {
        assertNull(afr.getCurrent());
        assertTrue(afr.next());

        assertNotNull(afr.getCurrent());
        final AlleleCounts ac = afr.getCurrent();
        assertEquals("A", ac.getReferenceAllele());
        assertEquals(5, ac.count("A"));
        assertEquals(2, ac.count("G"));
        assertEquals(0, ac.position());
        assertEquals("blah", afr.getCurrentReference());
      }
    }
  }

  public void testTextFormatWithSkipped() throws Exception {
    final MemoryPrintStream mps = new MemoryPrintStream();
    Diagnostic.setLogStream(mps.printStream());
    try (TestDirectory td = new TestDirectory()) {
      final String singleLine = "blah\t2\t.\t0\nblah\t5\tA\t5\n";

      final File out = new File(td, "snps.ac.gz");
      FileHelper.stringToGzFile(singleLine, out);

      new TabixIndexer(out, new File(td, "snps.ac.gz.tbi")).saveAlleleCountsIndex();

      try (AlleleCountsFileReader afr = AlleleCountsFileReader.openAlleleCountReader(out, null)) {
        assertNull(afr.getCurrent());
        assertTrue(afr.next());

        assertNotNull(afr.getCurrent());
        final AlleleCounts ac = afr.getCurrent();
        assertEquals("A", ac.getReferenceAllele());
        assertEquals(5, ac.count("A"));
        assertEquals(4, ac.position());
        assertEquals("blah", afr.getCurrentReference());
      }
    }
    Diagnostic.setLogStream();
  }

  private static final String HEADER = "##fileformat=VCFv4.1\n"
                              + "##INFO=<ID=AC,Number=.,Type=Integer,Description=\"Alternate Allele Count\">\n"
                              + "##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total Allele Count\">\n"
                              + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

  public void testVcfFormatACAN() throws Exception {
    try (TestDirectory td = new TestDirectory()) {
      final String singleLine = HEADER
        + "20\t60479\trs149529999\tC\t<DEL>\t100\tPASS\tAN=2184;AC=4,2\n" //this will be skipped, as it's an SV
        + "20\t60479\trs149529999\tC\tT,A\t100\tPASS\tAN=2184;AC=4,2\n";

      final File out = new File(td, "snps.vcf.gz");
      FileHelper.stringToGzFile(singleLine, out);

      new TabixIndexer(out, new File(td, "snps.vcf.gz.tbi")).saveVcfIndex();

      try (AlleleCountsFileReader afr = AlleleCountsFileReader.openAlleleCountReader(out, null)) {
        assertNull(afr.getCurrent());
        assertTrue(afr.next());

        assertNotNull(afr.getCurrent());
        final AlleleCounts ac = afr.getCurrent();
        assertEquals("C", ac.getReferenceAllele());
        assertEquals(4, ac.count("T"));
        assertEquals(2, ac.count("A"));
        assertEquals(-1, ac.count("G"));
        assertEquals(2178, ac.count("C"));
        assertEquals(60478, ac.position());
        assertEquals("20", afr.getCurrentReference());
      }
    }
  }
}
