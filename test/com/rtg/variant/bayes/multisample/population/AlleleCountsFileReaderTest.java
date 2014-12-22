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
