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

package com.rtg.visualization;

import static com.rtg.util.StringUtils.LS;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;

import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.io.FileUtils;
import com.rtg.util.test.FileHelper;
import com.rtg.vcf.header.VcfHeader;

import htsjdk.samtools.SAMRecord;

import junit.framework.TestCase;

/**
 */
public class AviewModelTest extends TestCase {

  private static final String MODEL_SNP = ""
    + VcfHeader.MINIMAL_HEADER + "\tSAMPLE" + LS
    + "g1" + "\t" + "1" + "\t" + "." + "\t" + "A" + "\t" + "C" + "\t" + "4.0" + "\t" + "PASS\t.\tGT\t0/1" + LS
    + "g1" + "\t" + "3" + "\t" + "." + "\t" + "G" + "\t" + "GAAA,GCC" + "\t" + "4.0" + "\t" + "PASS\t.\tGT\t1/2" + LS
    + "g1" + "\t" + "4" + "\t" + "." + "\t" + "AA" + "\t" + "AAA,AC" + "\t" + "4.0" + "\t" + "PASS\t.\tGT\t1/2" + LS
    ;

  private static final String MODEL_GENERATED_SNP = ""
    + VcfHeader.MINIMAL_HEADER + "\tSAMPLE" + LS
    + "g1" + "\t" + "1" + "\t" + "." + "\t" + "A" + "\t" + "C" + "\t" + "4.0" + "\t" + "PASS\t.\tGT\t0/1" + LS
    + "g1" + "\t" + "7" + "\t" + "." + "\t" + "G" + "\t" + "GCC,GAAA" + "\t" + "4.0" + "\t" + "PASS\t.\tGT\t1/2" + LS
    + "g1" + "\t" + "9" + "\t" + "." + "\t" + "AA" + "\t" + "AC,AAA" + "\t" + "4.0" + "\t" + "PASS\t.\tGT\t1/2" + LS
    ;

  /**
   * @throws IOException
   */
  public final void testReference() throws IOException {
    Diagnostic.setLogStream();
    final File f = FileUtils.createTempDir("aviewmodel", "test");
    try {

      AviewTest.prepareData(f, AviewTest.SAM);
      //System.out.println(AviewTest.SAM);
      //System.out.println(MODEL_SNP);
      //System.out.println(MODEL_GENERATED_SNP);
      final File snps = new File(f, "model.vcf");
      FileUtils.stringToFile(MODEL_SNP, snps);
      final File generatedsnps = new File(f, "model.generated.vcf");
      FileUtils.stringToFile(MODEL_GENERATED_SNP, generatedsnps);
      final AviewParams p = new AviewParamsBuilder()
        .alignments(new File[]{new File(f, AviewTest.ALIGNMENTS_BAM_FILE_NAME)})
        .trackFiles(new File[]{snps})
        .baselineFile(generatedsnps)
        .reference(new File(f, AviewTest.TEMPLATE_FILE_NAME))
        .region("g1:1-10")
        .create();

      final AviewModel m = new AviewModel(p);

      assertEquals(1, m.oneBasedStart());
      assertEquals(0, m.zeroBasedStart());
      assertEquals(13, m.zeroBasedEnd());
      final ArrayList<AviewVariant> called = ((AviewModel.CallSet) m.tracks().get(1)).sampleSnps().get("SAMPLE");
      final ArrayList<AviewVariant> generated = ((AviewModel.BaselineSet) m.tracks().get(0)).sampleSnps().get("SAMPLE");
      assertEquals(3, called.size());
      assertEquals(3, generated.size());
      final int[] inserts = m.inserts();
      assertEquals(13, inserts.length);
                                            //1, 2, 3, 4, 5, 6, 7, 8, 9, 0, 1, 2, 3
      final int[] expectedInserts = new int[]{0, 0, 0, 3, 1, 0, 0, 3, 0, 1, 0, 0, 0};
      //System.err.println(Arrays.toString(inserts));
      assertTrue(Arrays.equals(inserts, expectedInserts));

      assertEquals(5, m.records().size());

    } finally {
      AviewTest.deleteBrokenBam(f);
    }
  }

  /**
   * @throws IOException
   */
  public final void testSmallReference() throws IOException {
    Diagnostic.setLogStream();
    final File f = FileUtils.createTempDir("aviewmodel", "test");
    try {

      AviewTest.prepareData(f, AviewTest.SAM);
      //System.out.println(AviewTest.SAM);
      //System.out.println(MODEL_SNP);
      //System.out.println(MODEL_GENERATED_SNP);
      final File snps = new File(f, "model.vcf");
      FileUtils.stringToFile(MODEL_SNP, snps);
      final File generatedsnps = new File(f, "model.generated.vcf");
      FileUtils.stringToFile(MODEL_GENERATED_SNP, generatedsnps);
      final AviewParams p = new AviewParamsBuilder()
        .alignments(new File[]{new File(f, AviewTest.ALIGNMENTS_BAM_FILE_NAME)})
        .trackFiles(new File[]{snps})
        .baselineFile(generatedsnps)
        .reference(new File(f, AviewTest.TEMPLATE_FILE_NAME))
        .region("g1:4-5")
        .create();

      final AviewModel m = new AviewModel(p);

      assertEquals(3, m.oneBasedStart());
      assertEquals(2, m.zeroBasedStart());
      assertEquals(12, m.zeroBasedEnd());
      final ArrayList<AviewVariant> called = ((AviewModel.CallSet) m.tracks().get(1)).sampleSnps().get("SAMPLE");
      final ArrayList<AviewVariant> generated = ((AviewModel.BaselineSet) m.tracks().get(0)).sampleSnps().get("SAMPLE");
      assertEquals(2, called.size());
      assertEquals(2, generated.size());
      final int[] inserts = m.inserts();
      assertEquals(10, inserts.length);
                                            //3, 4, 5, 6, 7, 8, 9, 0, 1, 2
      final int[] expectedInserts = new int[]{0, 3, 1, 0, 0, 3, 0, 1, 0, 0};
      //System.err.println(Arrays.toString(inserts));
      assertTrue(Arrays.equals(inserts, expectedInserts));

      assertEquals(3, m.records().size());

    } finally {
      AviewTest.deleteBrokenBam(f);
    }
  }

  public final void testExceptionalCondition() throws IOException {
    Diagnostic.setLogStream();
    final File f = FileUtils.createTempDir("aviewmodel", "test");
    try {

      AviewTest.prepareData(f, AviewTest.SAM);
      final File snps = new File(f, "model.snps");
      FileUtils.stringToFile(MODEL_SNP, snps);

     try {
       new AviewParamsBuilder()
         .alignments(new File[]{new File(f, AviewTest.ALIGNMENTS_BAM_FILE_NAME)})
         .trackFiles(new File[]{snps})
         .reference(new File(f, AviewTest.TEMPLATE_FILE_NAME))
         .region("g1:25-10")
         .create();
       fail();
     } catch (final IllegalArgumentException ex) {
       // Expected
     }

     try {
       new AviewParamsBuilder()
         .alignments(new File[]{new File(f, AviewTest.ALIGNMENTS_BAM_FILE_NAME)})
         .trackFiles(new File[]{snps})
         .reference(new File(f, AviewTest.TEMPLATE_FILE_NAME))
         .region("g1:10-5")
         .create();
       fail();
     } catch (final IllegalArgumentException ex) {
       // Expected
     }
    } finally {
      AviewTest.deleteBrokenBam(f);
    }
  }

  public void testSort() throws Exception {
    final String t = "\t";
    final String sam = "" + "@HD" + t + "VN:1.0" + t + "SO:coordinate" + LS
    + "@SQ" + t + "SN:g1" + t + "LN:20" + LS
    + "@RG" + t + "ID:c" + t + "SM:b" + LS
    + "@RG" + t + "ID:d" + t + "SM:a" + LS
    + "4" + t + "0" + t + "g1" + t + "1" + t + "255" + t + "8M" + t + "*" + t + "0" + t + "0" + t + "CGACTGGT" + t + "````````" + t + "AS:i:1" + t + "IH:i:1" + t + "RG:Z:d" + LS
    + "0" + t + "0" + t + "g1" + t + "3" + t + "255" + t + "8M" + t + "*" + t + "0" + t + "0" + t + "ATCGACTG" + t + "&'(`````" + t + "AS:i:0" + t + "IH:i:1" + LS
    + "1" + t + "0" + t + "g1" + t + "3" + t + "255" + t + "8M" + t + "*" + t + "0" + t + "0" + t + "ATCGACTG" + t + "````````" + t + "AS:i:0" + t + "IH:i:1" + t + "RG:Z:d" + LS
    + "2" + t + "0" + t + "g1" + t + "5" + t + "255" + t + "8M" + t + "*" + t + "0" + t + "0" + t + "CGACTGGT" + t + "````````" + t + "AS:i:1" + t + "IH:i:1" + t + "RG:Z:c" + LS
    + "3" + t + "0" + t + "g1" + t + "7" + t + "255" + t + "8M" + t + "*" + t + "0" + t + "0" + t + "ATCGACTG" + t + "&'(`````" + t + "AS:i:0" + t + "IH:i:1" + LS;

    Diagnostic.setLogStream();
    final File f = FileUtils.createTempDir("aviewmodel", "test");
    try {

      AviewTest.prepareData(f, sam);
      final File snps = new File(f, "model.vcf");
      FileUtils.stringToFile(MODEL_SNP, snps);

      final AviewParams p2 = new AviewParamsBuilder()
      .alignments(new File[]{new File(f, AviewTest.ALIGNMENTS_BAM_FILE_NAME)})
      .trackFiles(new File[]{snps})
      .sortReadGroup(true)
      .reference(new File(f, AviewTest.TEMPLATE_FILE_NAME))
      .region("g1:1-10")
      .create();

      final AviewModel avm = new AviewModel(p2);

      final ArrayList<SAMRecord> recs = avm.records();
      assertEquals(5, recs.size());
      assertEquals("2", recs.get(0).getReadName());
      assertEquals("4", recs.get(1).getReadName());
      assertEquals("1", recs.get(2).getReadName());
      assertEquals("0", recs.get(3).getReadName());
      assertEquals("3", recs.get(4).getReadName());

    } finally {
      AviewTest.deleteBrokenBam(f);
    }
  }


  public void testInsertsAtEndOfRegion() throws Exception {
    final File tmpDir = FileHelper.createTempDirectory();
    try {
      final String t = "\t";
      final String sam = "" + "@HD" + t + "VN:1.0" + t + "SO:coordinate" + LS
          + "@SQ" + t + "SN:g1" + t + "LN:20" + LS
          + "@RG" + t + "ID:c" + t + "SM:b" + LS
          + "@RG" + t + "ID:d" + t + "SM:a" + LS
          + "0" + t + "0" + t + "g1" + t + "0" + t + "255" + t + "8M" + t + "*" + t + "0" + t + "0" + t + "AAATCGAC" + t + "````````" + t + "AS:i:0" + t + "IH:i:1" + t + "RG:Z:d" + LS
          + "1" + t + "0" + t + "g1" + t + "1" + t + "255" + t + "7M2I1M" + t + "*" + t + "0" + t + "0" + t + "AATCGACAAT" + t + "&'(```````" + t + "AS:i:3" + t + "IH:i:1" + LS;
      AviewTest.prepareData(tmpDir, sam);

      final AviewParams p2 = new AviewParamsBuilder()
      .alignments(new File[]{new File(tmpDir, AviewTest.ALIGNMENTS_BAM_FILE_NAME)})
      .sortReadGroup(true)
      .reference(new File(tmpDir, AviewTest.TEMPLATE_FILE_NAME))
      .region("g1:5-8")
      .create();

      final AviewModel avm = new AviewModel(p2);  //this merely passing tests the expected behaviour.

      assertEquals(2, avm.inserts()[7]);
    } finally {
      FileHelper.deleteAll(tmpDir);
    }
  }
}
