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
package com.rtg.relation;

import java.io.File;
import java.io.IOException;

import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.io.MemoryPrintStream;
import com.rtg.util.io.TestDirectory;
import com.rtg.util.test.FileHelper;
import com.rtg.util.test.NanoRegression;
import com.rtg.vcf.header.PedigreeField;
import com.rtg.vcf.header.VcfHeader;

import junit.framework.TestCase;

/**
 */
public class VcfPedigreeParserTest extends TestCase {

  private NanoRegression mNano = null;

  @Override
  public void setUp() {
    mNano = new NanoRegression(this.getClass());
  }
  @Override
  public void tearDown() throws Exception {
    try {
      mNano.finish();
    } finally {
      mNano = null;
    }
  }

  public void testParsing() throws IOException {
    try (final TestDirectory dir = new TestDirectory()) {
      final File vcfFile = FileHelper.resourceToFile("com/rtg/relation/resources/vcfheader.vcf", new File(dir, "header.vcf"));
      final GenomeRelationships ped2 = VcfPedigreeParser.loadFile(vcfFile);
      mNano.check("pedfromvcf", PedFileParser.toString(ped2));
    }
  }

  public void testPedConversion() throws IOException {
    try (final TestDirectory dir = new TestDirectory()) {
      final File pedFile = FileHelper.resourceToFile("com/rtg/relation/resources/pop.ped", new File(dir, "pop.ped"));
      final GenomeRelationships ped = PedFileParser.loadFile(pedFile);
      final VcfHeader header = new VcfHeader();
      header.addLine(VcfHeader.VERSION_LINE);
      VcfPedigreeParser.addPedigreeFields(header, ped);
      for (String sample : ped.filter(new GenomeRelationships.PrimaryGenomeFilter(ped)).genomes()) {
        header.addSampleName(sample);
      }
      mNano.check("vcffromped.vcf", header.toString());
    }
  }

  public void testFullCircle() throws IOException {
    try (final TestDirectory dir = new TestDirectory()) {
      final File vcfFile = FileHelper.resourceToFile("com/rtg/relation/resources/vcffromped.vcf", new File(dir, "pop.vcf"));
      final GenomeRelationships ped = VcfPedigreeParser.loadFile(vcfFile);
      final VcfHeader header = new VcfHeader();
      header.addLine(VcfHeader.VERSION_LINE);
      VcfPedigreeParser.addPedigreeFields(header, ped);
      for (String sample : ped.filter(new GenomeRelationships.PrimaryGenomeFilter(ped)).genomes()) {
        header.addSampleName(sample);
      }
      mNano.check("vcffromped.vcf", header.toString());
    }
  }

  public void testFullCircleWithCellLine() throws IOException {
    try (final TestDirectory dir = new TestDirectory()) {
      final File vcfFile = FileHelper.resourceToFile("com/rtg/relation/resources/derived.vcf", new File(dir, "pop.vcf"));
      final GenomeRelationships ped = VcfPedigreeParser.loadFile(vcfFile);
      final VcfHeader header = new VcfHeader();
      header.addLine(VcfHeader.VERSION_LINE);
      VcfPedigreeParser.addPedigreeFields(header, ped);
      for (String sample : ped.filter(new GenomeRelationships.PrimaryGenomeFilter(ped)).genomes()) {
        header.addSampleName(sample);
      }
      mNano.check("derived.vcf", header.toString());
    }
  }

  public void testNoInformationWarning() throws IOException {
    final GenomeRelationships ped = new GenomeRelationships();
    final PedigreeField f = new PedigreeField("##PEDIGREE=<Sibling=>");
    MemoryPrintStream mps = new MemoryPrintStream();
    Diagnostic.setLogStream(mps.printStream());
    VcfPedigreeParser.parsePedLine(ped, f);
    assertTrue(mps.toString().contains("Pedigree line contains no pedigree information: ##PEDIGREE=<Sibling=>"));
    Diagnostic.setLogStream();
  }

}
