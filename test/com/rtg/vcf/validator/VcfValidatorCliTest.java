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

package com.rtg.vcf.validator;

import java.io.File;
import java.io.IOException;

import com.rtg.launcher.AbstractCli;
import com.rtg.launcher.AbstractCliTest;
import com.rtg.util.io.TestDirectory;
import com.rtg.util.test.FileHelper;

/**
 */
public class VcfValidatorCliTest extends AbstractCliTest {

  @Override
  protected AbstractCli getCli() {
    return new VcfValidatorCli();
  }

  public void testHelp() {
    checkHelp("rtg vcfvalidator [OPTION]... FILE"
        , "Validates the contents of a VCF file conform to expected value ranges."
        , "VCF format file to be validated"
        );
    checkExtendedHelp("rtg vcfvalidator [OPTION]... FILE"
        , "--Xrules=FILE", "File defining rules for validation of VCF input"
        , "--Xverbose", "Set to output all failed records to error output instead of the first 10"
        );
  }

  public void testNormalNano() throws IOException {
    try (TestDirectory dir = new TestDirectory()) {
      final File vcf = new File(dir, "fullOfBad.vcf");
      FileHelper.resourceToFile("com/rtg/vcf/validator/resources/fullOfBad.vcf", vcf);
      final String errOut = checkMainInitBadFlags(vcf.getPath());
      mNano.check("fullOfBad.txt", errOut, true);
    }
  }

  public void testVerboseNano() throws IOException {
    try (TestDirectory dir = new TestDirectory()) {
      final File vcf = new File(dir, "fullOfBad.vcf");
      FileHelper.resourceToFile("com/rtg/vcf/validator/resources/fullOfBad.vcf", vcf);
      final String errOut = checkMainInitBadFlags("--Xverbose", vcf.getPath());
      mNano.check("fullOfBad_verbose.txt", errOut, true);
    }
  }

  public void testMismatchedHeaders() throws IOException {
    try (TestDirectory dir = new TestDirectory()) {
      final File vcf = new File(dir, "badHeaders.vcf");
      FileHelper.resourceToFile("com/rtg/vcf/validator/resources/badHeaders.vcf", vcf);
      final String errOut = checkMainInitBadFlags(vcf.getPath());
      mNano.check("badHeaders.txt", errOut, true);
    }
  }

  public void testBoringNormal() throws IOException {
    try (TestDirectory dir = new TestDirectory()) {
      final File vcf = new File(dir, "boringNormal.vcf");
      FileHelper.resourceToFile("com/rtg/vcf/validator/resources/boringNormal.vcf", vcf);
      final String out = checkMainInitOk(vcf.getPath());
      assertEquals("", out);
    }
  }

}
