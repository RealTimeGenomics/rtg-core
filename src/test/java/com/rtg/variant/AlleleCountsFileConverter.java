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
