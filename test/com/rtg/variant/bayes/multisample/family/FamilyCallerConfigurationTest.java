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

package com.rtg.variant.bayes.multisample.family;

import java.io.BufferedReader;
import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.StringReader;

import com.rtg.reference.Ploidy;
import com.rtg.relation.GenomeRelationships;
import com.rtg.util.StringUtils;
import com.rtg.util.TestUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.NoTalkbackSlimException;
import com.rtg.util.test.BgzipFileHelper;
import com.rtg.util.test.FileHelper;
import com.rtg.variant.AlleleCountsFileConverter;
import com.rtg.variant.GenomePriorParams;
import com.rtg.variant.VariantParams;
import com.rtg.variant.VariantParamsBuilder;
import com.rtg.variant.format.VariantOutputVcfFormatter;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMReadGroupRecord;

import junit.framework.TestCase;

/**
 */
public class FamilyCallerConfigurationTest extends TestCase {

  static final String TRIO_PED = ""
      + "0\tfather\t0\t0\t1\t0" + StringUtils.LS
      + "0\tmother\t0\t0\t2\t0" + StringUtils.LS
      + "0\tchild\tfather\tmother\t1\t0" + StringUtils.LS;

  static final String TOOMANY_PED = TRIO_PED
      + "0\tchild2\tfather\tanothermother\t1\t0" + StringUtils.LS;


  public void testCreation() throws Exception {
    final File tmp = FileHelper.createTempDirectory();
    try {
      final File popFile = new File(tmp, "pop.vcf.gz");
      BgzipFileHelper.resourceToBgzipFile("com/rtg/variant/resources/pop58.vcf", popFile);

      final File alleleCountFile = new File(tmp, "allele.ac");

      final AlleleCountsFileConverter blah = new AlleleCountsFileConverter();
      blah.convert(popFile, alleleCountFile);

      Diagnostic.setLogStream();

      try {
        final GenomeRelationships rel = GenomeRelationships.loadGenomeRelationships(new BufferedReader(new StringReader(TOOMANY_PED)));
        final VariantParams p = new VariantParamsBuilder().genomeRelationships(rel).create();
        new FamilyCallerConfiguration.Configurator().getConfig(p, new String[] {"child", "mother", "father"});
        fail("Should not have liked the family");
      } catch (NoTalkbackSlimException e) {
        // Expected
      }

      final GenomeRelationships rel = GenomeRelationships.loadGenomeRelationships(new BufferedReader(new StringReader(TRIO_PED)));
      final VariantParamsBuilder b = new VariantParamsBuilder();
      b.genomePriors(GenomePriorParams.builder().create());
      b.genomeRelationships(rel);
      b.machineErrorName("illumina");
      b.populationPriors(alleleCountFile);
      final VariantParams p = b.create();

      final FamilyCallerConfiguration config = new FamilyCallerConfiguration.Configurator().getConfig(p, new String[] {"child", "mother", "father"});
      assertNotNull(config.getGenomeNames());
      assertNotNull(config.getJointCaller());
      assertEquals(3, config.numberOfGenomes());
      assertEquals(3, config.getGenomeNames().length);
      assertEquals("father", config.getGenomeNames()[0]);
      assertEquals("mother", config.getGenomeNames()[1]);
      assertEquals("child", config.getGenomeNames()[2]);
      assertTrue(config.getJointCaller() instanceof FamilyCaller);
      assertTrue(config.handlesPloidy(Ploidy.POLYPLOID));
      assertTrue(config.handlesPloidy(Ploidy.HAPLOID));
      assertTrue(config.handlesPloidy(Ploidy.DIPLOID));

      final VariantOutputVcfFormatter output = config.getOutputFormatter(p);
      final ByteArrayOutputStream bos = new ByteArrayOutputStream();
      final SAMFileHeader sfh = new SAMFileHeader();
      final SAMReadGroupRecord samReadGroupRecord = new SAMReadGroupRecord("ab");
      samReadGroupRecord.setSample("father");
      final SAMReadGroupRecord samReadGroupRecord2 = new SAMReadGroupRecord("abc");
      samReadGroupRecord2.setSample("mother");
      final SAMReadGroupRecord samReadGroupRecord3 = new SAMReadGroupRecord("abcd");
      samReadGroupRecord3.setSample("child");
      sfh.addReadGroup(samReadGroupRecord);
      sfh.addReadGroup(samReadGroupRecord2);
      sfh.addReadGroup(samReadGroupRecord3);
      output.writeHeader(bos, p, sfh);
      final String header = bos.toString();
      //System.err.println(header);
      TestUtils.containsAll(header, "##PEDIGREE=<Child=child,Mother=mother,Father=father>");
    } finally {
      FileHelper.deleteAll(tmp);
    }
  }

  public void testBadRelationships() throws Exception {
    Diagnostic.setLogStream();

    final GenomeRelationships rel = GenomeRelationships.loadGenomeRelationships(new BufferedReader(new StringReader(TRIO_PED)));
    final VariantParamsBuilder b = new VariantParamsBuilder();
    b.genomePriors(GenomePriorParams.builder().create());
    b.genomeRelationships(rel);
    b.machineErrorName("illumina");
    final VariantParams p = b.create();

    FamilyCallerConfiguration config = new FamilyCallerConfiguration.Configurator().getConfig(p, new String[] {"child", "mother", "father"});
    assertEquals(3, config.numberOfGenomes());

    config = new FamilyCallerConfiguration.Configurator().getConfig(p, new String[] {"child", "mother"});
    assertEquals(3, config.numberOfGenomes());

    config = new FamilyCallerConfiguration.Configurator().getConfig(p, new String[] {"child", "father"});
    assertEquals(3, config.numberOfGenomes());

    try {
      new FamilyCallerConfiguration.Configurator().getConfig(p, new String[] {"father", "mother"});
      fail("Accepted childless family");
    } catch (NoTalkbackSlimException e) {
      assertEquals("Not enough family members have mapping data provided", e.getMessage());
    }

    try {
      new FamilyCallerConfiguration.Configurator().getConfig(p, new String[] {"father"});
      fail("Accepted father only family");
    } catch (NoTalkbackSlimException e) {
      assertEquals("Not enough family members have mapping data provided", e.getMessage());
    }

    try {
      new FamilyCallerConfiguration.Configurator().getConfig(p, new String[] {"child"});
      fail("Accepted child only family");
    } catch (NoTalkbackSlimException e) {
      assertEquals("Not enough family members have mapping data provided", e.getMessage());
    }
  }
}
