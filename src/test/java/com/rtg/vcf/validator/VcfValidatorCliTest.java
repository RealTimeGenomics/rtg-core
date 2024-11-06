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
import com.rtg.variant.format.VcfInfoField;
import com.rtg.vcf.DefaultVcfWriter;
import com.rtg.vcf.VcfRecord;
import com.rtg.vcf.VcfWriter;
import com.rtg.vcf.header.FormatField;
import com.rtg.vcf.header.InfoField;
import com.rtg.vcf.header.MetaType;
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
   */
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
          final String val = getTypeValue(formatField.getType());
          vcfRecord.addFormatAndSample(formatField.getId(), val);
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

  /**
   * check that all declared INFO fields are present in the rules file
   * should catch missing rules earlier than regression tests
   */
  public void testAllInfoFieldsHaveARule() throws IOException {
    try (TestDirectory td = new TestDirectory()) {
      final VcfHeader header = new VcfHeader();
      header.setVersionValue(VcfHeader.VERSION_VALUE);
      for (VcfInfoField vcfInfoField : VcfInfoField.values()) {
        vcfInfoField.updateHeader(header);
      }
      final File file = new File(td, "foo.vcf");
      try (final OutputStream out = new FileOutputStream(file)) {
        final VcfWriter writer = new DefaultVcfWriter(header, out);
        final VcfRecord vcfRecord = new VcfRecord("foo", 3, "A");
        vcfRecord.setNumberOfSamples(0);
        for (InfoField infoField : header.getInfoLines()) {
          final String val = getTypeValue(infoField.getType());
          vcfRecord.setInfo(infoField.getId(), val);
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

  private String getTypeValue(MetaType metaType) {
    switch (metaType) {
      case CHARACTER:
        return ".";
      case STRING:
        return "AB";
      case INTEGER:
        return "1";
      case FLOAT:
        return "-0.11";
      case FLAG:
        return ".";
      default:
        return ".";
    }
  }

}
