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
package com.rtg.vcf;

import java.io.File;

import com.rtg.launcher.AbstractCli;
import com.rtg.launcher.AbstractCliTest;
import com.rtg.util.StringUtils;
import com.rtg.util.TestUtils;
import com.rtg.util.io.TestDirectory;
import com.rtg.util.test.FileHelper;

/**
 */
public class VcfSubsetTest extends AbstractCliTest {

  public void testFlags() {
    checkHelp("rtg vcfsubset", "Removes columnar data from VCF records");
  }

  @Override
  protected AbstractCli getCli() {
    return new VcfSubset();
  }

  public void testKeepInfoACAN() throws Exception {
    try (TestDirectory td = new TestDirectory()) {
      final File f = FileHelper.resourceToGzFile("com/rtg/vcf/resources/vcfsubset.vcf", new File(td, "vcf.vcf.gz"));
      final File out = new File(td, "out.vcf");

      checkMainInitOk("-i", f.getPath(), "-o", out.getPath(), "--keep-info", "AC", "--keep-info", "AN", "-Z");

      final String content = FileHelper.fileToString(out);
      final String editedContent = StringUtils.grepMinusV(StringUtils.grepMinusV(content, "^##RUN-ID"), "^##CL");
      mNano.check("vcfsubset-keepinfoACAN.vcf", editedContent);
    }
  }

  public void testKeepFilter() throws Exception {
    try (TestDirectory td = new TestDirectory()) {
      final File f = FileHelper.resourceToGzFile("com/rtg/vcf/resources/vcfsubset.vcf", new File(td, "vcf.vcf.gz"));
      final File out = new File(td, "out.vcf");

      checkMainInitOk("-i", f.getPath(), "-o", out.getPath(), "--keep-filter", "YEA", "-Z");

      final String content = FileHelper.fileToString(out);
      final String editedContent = StringUtils.grepMinusV(StringUtils.grepMinusV(content, "^##RUN-ID"), "^##CL");
      mNano.check("vcfsubset-keepfilter.vcf", editedContent);
    }
  }

  public void testKeepSamples() throws Exception {
    try (TestDirectory td = new TestDirectory()) {
      final File f = FileHelper.resourceToGzFile("com/rtg/vcf/resources/vcfsubset.vcf", new File(td, "vcf.vcf.gz"));
      final File out = new File(td, "out.vcf");

      checkMainInitOk("-i", f.getPath(), "-o", out.getPath(), "--keep-sample", "HG00096", "--keep-sample", "HG00100", "-Z");

      final String content = FileHelper.fileToString(out);
      final String editedContent = StringUtils.grepMinusV(StringUtils.grepMinusV(content, "^##RUN-ID"), "^##CL");
      mNano.check("vcfsubset-keepsamples.vcf", editedContent);
    }
  }

  public void testRemoveFormat() throws Exception {
    try (TestDirectory td = new TestDirectory()) {
      final File f = FileHelper.resourceToGzFile("com/rtg/vcf/resources/vcfsubset.vcf", new File(td, "vcf.vcf.gz"));
      final File out = new File(td, "out.vcf");

      checkMainInitOk("-i", f.getPath(), "-o", out.getPath(), "--remove-format", "DS", "-Z");

      final String content = FileHelper.fileToString(out);
      final String editedContent = StringUtils.grepMinusV(StringUtils.grepMinusV(content, "^##RUN-ID"), "^##CL");
      mNano.check("vcfsubset-removeformat.vcf", editedContent);
    }
  }

  public void testRemoveMulti() throws Exception {
    try (TestDirectory td = new TestDirectory()) {
      final File f = FileHelper.resourceToGzFile("com/rtg/vcf/resources/vcfsubset.vcf", new File(td, "vcf.vcf.gz"));
      final File out = new File(td, "out.vcf");

      checkMainInitOk("-i", f.getPath(), "-o", out.getPath(), "--remove-samples", "--keep-info", "AN", "--keep-info", "AC", "--keep-filter", "YEA", "-Z");

      final String content = FileHelper.fileToString(out);
      final String editedContent = StringUtils.grepMinusV(StringUtils.grepMinusV(content, "^##RUN-ID"), "^##CL");
      mNano.check("vcfsubset-multi.vcf", editedContent);
    }
  }

  public void testValidation() throws Exception {
    try (TestDirectory main = new TestDirectory()) {
      final File in = FileHelper.resourceToGzFile("com/rtg/vcf/resources/vcfsubset.vcf", new File(main, "vcf.vcf.gz"));
      final File out = new File(main, "out.gz");
      final File out2 = new File(main, "out.txt");
      String err = checkHandleFlagsErr("-i", "foo", "-o", out.getPath());
      TestUtils.containsAll(err, "Given file \"foo\" does not exist.");

      err = checkHandleFlagsErr("-i", main.getPath(), "-o", out.getPath());
      TestUtils.containsAll(TestUtils.unwrap(err), main.getPath() + "\" is a directory");

      assertTrue(out.createNewFile());
      err = checkHandleFlagsErr("-i", in.getPath(), "-o", new File(main, "out").getPath());
      TestUtils.containsAll(TestUtils.unwrap(err), "The file \"" + out.getPath() + "\" already exists");

      assertTrue(out2.createNewFile());
      err = checkHandleFlagsErr("-i", in.getPath(), "-o", out2.getPath(), "--no-gzip");
      TestUtils.containsAll(TestUtils.unwrap(err), "The file \"" + out2.getPath() + "\" already exists");

      err = checkHandleFlagsErr("-i", in.getPath(), "-o", new File(main, "newout.gz").getPath(), "--remove-infos", "--remove-info", "blah");
      TestUtils.containsAll(TestUtils.unwrap(err), "Only one of --remove-infos, --remove-info, or --keep-info can be set");
      err = checkHandleFlagsErr("-i", in.getPath(), "-o", new File(main, "newout.gz").getPath(), "--remove-infos", "--keep-info", "blah");
      TestUtils.containsAll(TestUtils.unwrap(err), "Only one of --remove-infos, --remove-info, or --keep-info can be set");
      err = checkHandleFlagsErr("-i", in.getPath(), "-o", new File(main, "newout.gz").getPath(), "--remove-info", "feh", "--keep-info", "blah");
      TestUtils.containsAll(TestUtils.unwrap(err), "Only one of --remove-infos, --remove-info, or --keep-info can be set");

      err = checkHandleFlagsErr("-i", in.getPath(), "-o", new File(main, "newout.gz").getPath(), "--remove-filters", "--remove-filter", "blah");
      TestUtils.containsAll(TestUtils.unwrap(err), "Only one of --remove-filters, --remove-filter, or --keep-filter can be set");
      err = checkHandleFlagsErr("-i", in.getPath(), "-o", new File(main, "newout.gz").getPath(), "--remove-filters", "--keep-filter", "blah");
      TestUtils.containsAll(TestUtils.unwrap(err), "Only one of --remove-filters, --remove-filter, or --keep-filter can be set");
      err = checkHandleFlagsErr("-i", in.getPath(), "-o", new File(main, "newout.gz").getPath(), "--remove-filter", "feh", "--keep-filter", "blah");
      TestUtils.containsAll(TestUtils.unwrap(err), "Only one of --remove-filters, --remove-filter, or --keep-filter can be set");

      err = checkHandleFlagsErr("-i", in.getPath(), "-o", new File(main, "newout.gz").getPath(), "--remove-samples", "--remove-sample", "blah");
      TestUtils.containsAll(TestUtils.unwrap(err), "Only one of --remove-samples, --remove-sample, or --keep-sample can be set");
      err = checkHandleFlagsErr("-i", in.getPath(), "-o", new File(main, "newout.gz").getPath(), "--remove-samples", "--keep-sample", "blah");
      TestUtils.containsAll(TestUtils.unwrap(err), "Only one of --remove-samples, --remove-sample, or --keep-sample can be set");
      err = checkHandleFlagsErr("-i", in.getPath(), "-o", new File(main, "newout.gz").getPath(), "--remove-sample", "feh", "--keep-sample", "blah");
      TestUtils.containsAll(TestUtils.unwrap(err), "Only one of --remove-samples, --remove-sample, or --keep-sample can be set");

      err = checkHandleFlagsErr("-i", in.getPath(), "-o", new File(main, "newout.gz").getPath(), "--remove-format", "feh", "--keep-format", "blah");
      TestUtils.containsAll(TestUtils.unwrap(err), "Only one of --remove-format or --keep-format can be set");
    }
  }

  public void testMissingSample() throws Exception {
    try (TestDirectory td = new TestDirectory()) {
      final File f = FileHelper.resourceToGzFile("com/rtg/vcf/resources/vcfsubset.vcf", new File(td, "vcf.vcf.gz"));
      final File out = new File(td, "out.vcf");

      assertEquals("Error: Sample fields not contained in VCF header: BL RJ" + StringUtils.LS, checkMainInitBadFlags("-i", f.getPath(), "-o", out.getPath(), "--keep-sample", "HG00097", "--keep-sample", "HG00099", "--keep-sample", "BL", "--keep-sample", "RJ", "-Z"));
    }
  }

  public void testMissingInfo() throws Exception {
    try (TestDirectory td = new TestDirectory()) {
      final File f = FileHelper.resourceToGzFile("com/rtg/vcf/resources/vcfsubset.vcf", new File(td, "vcf.vcf.gz"));
      final File out = new File(td, "out.vcf");

      assertEquals("Error: Info fields not contained in VCF meta-information: BL RJ" + StringUtils.LS, checkMainInitBadFlags("-i", f.getPath(), "-o", out.getPath(), "--keep-info", "BL", "--keep-info", "RJ", "-Z"));
    }
  }

  public void testMissingFilter() throws Exception {
    try (TestDirectory td = new TestDirectory()) {
      final File f = FileHelper.resourceToGzFile("com/rtg/vcf/resources/vcfsubset.vcf", new File(td, "vcf.vcf.gz"));
      final File out = new File(td, "out.vcf");

      assertEquals("Error: Filter fields not contained in VCF meta-information: BL RJ" + StringUtils.LS, checkMainInitBadFlags("-i", f.getPath(), "-o", out.getPath(), "--keep-filter", "BL", "--keep-filter", "RJ", "-Z"));
    }
  }

  public void testMissingFormat() throws Exception {
    try (TestDirectory td = new TestDirectory()) {
      final File f = FileHelper.resourceToGzFile("com/rtg/vcf/resources/vcfsubset.vcf", new File(td, "vcf.vcf.gz"));
      final File out = new File(td, "out.vcf");

      assertEquals("Error: Format fields not contained in VCF meta-information: BL RJ" + StringUtils.LS, checkMainInitBadFlags("-i", f.getPath(), "-o", out.getPath(), "--keep-format", "BL", "--keep-format", "RJ", "-Z"));
    }
  }

  public void testExplosion() throws Exception {
    try (TestDirectory td = new TestDirectory()) {
      final File f = FileHelper.resourceToGzFile("com/rtg/vcf/resources/vcfsubset.vcf", new File(td, "vcf.vcf.gz"));
      final File out = new File(td, "out.vcf");

      assertEquals("Records skipped due to invalid sample fields: 2" + StringUtils.LS, checkMainInitWarn("-i", f.getPath(), "-o", out.getPath(), "--remove-format", "GT", "--remove-format", "DS", "--remove-format", "GL", "-Z"));

      final String content = FileHelper.fileToString(out);
      assertEquals("", StringUtils.grepMinusV(content, "^#")); //all the records get wiped out by keeping
    }
  }
}
