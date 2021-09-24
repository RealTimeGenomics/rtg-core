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
package com.rtg.variant;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;

import com.reeltwo.jumble.annotations.TestClass;
import com.rtg.tabix.TabixIndexer;
import com.rtg.tabix.UnindexableDataException;
import com.rtg.util.StringUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.io.FileUtils;
import com.rtg.variant.bayes.multisample.population.AlleleCounts;
import com.rtg.variant.bayes.multisample.population.AlleleCountsFileReader;


/**
 * Converts a <code>VCF</code> file containing a population into a
 * tab separated population allele counts file of the format:
 *
 * reference_name reference_pos [allele allele_count]+
 *
 * Where the first allele is that of the reference.
 */
@TestClass("com.rtg.variant.bayes.multisample.population.PopulationHwHypothesesCreatorTest")
public class AlleleCountsFileConverter {

  /**
   * Converts a <code>VCF</code> file containing a population into an allele counts file
   * @param vcfFile the <code>VCF</code> file to convert
   * @param out the file to output to
   * @throws IOException if something bad happen
   */
  public void convert(File vcfFile, File out) throws IOException {
    Diagnostic.developerLog("Loading variants...");
    try (OutputStream os = FileUtils.createOutputStream(out, true)) {
      try (AlleleCountsFileReader acr = AlleleCountsFileReader.openAlleleCountReader(vcfFile, null)) {
        while (acr.next()) {
          final AlleleCounts alleleCounts = acr.getCurrent();
          os.write(acr.getCurrentReference().getBytes());
          os.write('\t');
          os.write(Integer.toString(alleleCounts.position() + 1).getBytes());
          os.write('\t');

          final String refAllele = alleleCounts.getReferenceAllele();
          os.write(refAllele.getBytes());
          os.write('\t');
          os.write(Integer.toString(alleleCounts.count(refAllele)).getBytes());

          for (String allele : alleleCounts.allelesSeen()) {
            if (!allele.equals(refAllele)) {
              os.write('\t');
              os.write(allele.getBytes());
              os.write('\t');
              os.write(Integer.toString(alleleCounts.count(allele)).getBytes());
            }
          }
          os.write(StringUtils.LS.getBytes());
        }
      }
    }
    try {
      new TabixIndexer(out).saveAlleleCountsIndex();
    } catch (final UnindexableDataException e) {
      System.err.println("Cannot produce TABIX index for: " + out + ": " + e.getMessage());
    }
  }

  /**
   * @param args command line arguments
   * @throws Exception if rargh blargh
   */
  public static void main(String[] args) throws Exception {
    if (args.length != 2) {
      System.err.println("Supply [input vcf file] [output file] if you wish this to do something");
      System.exit(1);
    }
    Diagnostic.setLogStream(System.out);
    final AlleleCountsFileConverter ppc = new AlleleCountsFileConverter();
    ppc.convert(new File(args[0]), new File(args[1]));
  }
}
