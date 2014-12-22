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
package com.rtg.variant.util;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;

import com.rtg.util.Resources;
import com.rtg.util.io.MemoryPrintStream;
import com.rtg.util.test.NanoRegression;
import com.rtg.vcf.VcfReader;
import com.rtg.vcf.VcfRecord;
import com.rtg.vcf.VcfUtils;

import junit.framework.TestCase;

/**
 */
public class GenotypeProportionsTest extends TestCase {

  private NanoRegression mNano = null;
  @Override
  public void setUp() throws Exception {
    super.setUp();
    mNano = new NanoRegression(GenotypeProportionsTest.class);
  }

  @Override
  public void tearDown() throws Exception {
    super.tearDown();
    try {
      mNano.finish();
    } finally {
      mNano = null;
    }
  }

  public void test() throws IOException {
    final String vcfResource = "com/rtg/variant/util/resources/gt_prop_test.vcf";
    final String resultName = "gt_prop_results.txt";
    check(vcfResource, resultName);
  }

  public void test2() throws IOException {
    final String vcfResource = "com/rtg/variant/util/resources/gt_prop_ploidy_test.vcf";
    final String resultName = "gt_prop_ploidy_results.txt";
    check(vcfResource, resultName);
  }

  private void check(String vcfResource, String resultName) throws IOException {
    final GenotypeProportions prop = new GenotypeProportions();
    try (VcfReader r = new VcfReader(new BufferedReader(new InputStreamReader(Resources.getResourceAsStream(vcfResource))))) {
      while (r.hasNext()) {
        final VcfRecord rec = r.next();
        final ArrayList<String> sampleGts = rec.getFormatAndSample().get(VcfUtils.FORMAT_GENOTYPE);
        prop.addRecord(new Genotype(sampleGts.get(0)), new Genotype(sampleGts.get(1)), new Genotype(sampleGts.get(2)));
      }
    }
    MemoryPrintStream mps = new MemoryPrintStream();
    prop.writeResults(mps.printStream());
    mNano.check(resultName, mps.toString());
  }
}
