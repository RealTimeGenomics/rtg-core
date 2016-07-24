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

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;

import com.rtg.launcher.AbstractCli;
import com.rtg.launcher.AbstractCliTest;
import com.rtg.util.StringUtils;
import com.rtg.util.io.TestDirectory;
import com.rtg.util.test.FileHelper;
import com.rtg.variant.format.VcfFormatField;
import com.rtg.vcf.DefaultVcfWriter;
import com.rtg.vcf.VcfRecord;
import com.rtg.vcf.VcfWriter;
import com.rtg.vcf.header.FormatField;
import com.rtg.vcf.header.VcfHeader;

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

  /**
   * check that all declared format fields are present in the rules file
   * should catch missing rules earlier than regression tests
   * */
  public void testAllFormatFieldsHaveARule() throws IOException {
    try (TestDirectory td = new TestDirectory()) {
      final VcfHeader header = new VcfHeader();
      header.setVersionValue(VcfHeader.VERSION_VALUE);
      for (VcfFormatField vcfFormatField : VcfFormatField.values()) {
        vcfFormatField.updateHeader(header);
      }
      header.addSampleName("bar");
      final File file = new File(td, "foo.vcf");
      try (final OutputStream out = new FileOutputStream(file)) {
        final VcfWriter writer = new DefaultVcfWriter(header, out);
        final VcfRecord vcfRecord = new VcfRecord("foo", 3, "A");
        vcfRecord.setNumberOfSamples(1);
        for (FormatField formatField : header.getFormatLines()) {
          final VcfFormatField vcfFormatField = VcfFormatField.valueOf(formatField.getId());
          String val;
          switch (formatField.getType()) {
            case CHARACTER:
              val = ".";
              break;
            case STRING:
              val = "AB";
              break;
            case INTEGER:
              val = "1";
              break;
            case FLOAT:
              val = "-0.11";
              break;
            case FLAG:
              val = ".";
              break;
            default:
              val = ".";
          }
          vcfRecord.addFormatAndSample(vcfFormatField.name(), val);
        }
        writer.write(vcfRecord);
        writer.close();
      }
      final VcfValidatorCli cli = new VcfValidatorCli();
      final ByteArrayOutputStream out = new ByteArrayOutputStream();
      final ByteArrayOutputStream err = new ByteArrayOutputStream();
      try (PrintStream printErr = new PrintStream(err)) {
        cli.mainInit(new String[]{file.getAbsolutePath()}, out, printErr);
      }
      for (String line : StringUtils.split(err.toString(), '\n')) {
        assertFalse(line, line.contains("not contained in the specified rule set."));
      }
    }
  }

}
