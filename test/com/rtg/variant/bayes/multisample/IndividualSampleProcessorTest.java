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
package com.rtg.variant.bayes.multisample;

import static com.rtg.sam.SharedSamConstants.OUT_SAM;
import static com.rtg.sam.SharedSamConstants.REF_SEQS;
import static com.rtg.sam.SharedSamConstants.SAM1;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;

import com.rtg.launcher.OutputParams;
import com.rtg.launcher.SequenceParams;
import com.rtg.mode.SequenceMode;
import com.rtg.reader.ReaderTestUtils;
import com.rtg.reference.Sex;
import com.rtg.util.InvalidParamsException;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.io.FileUtils;
import com.rtg.util.test.FileHelper;
import com.rtg.variant.DefaultMachineErrorChooser;
import com.rtg.variant.GenomePriorParams;
import com.rtg.reference.SexMemo;
import com.rtg.variant.VariantAlignmentRecord;
import com.rtg.variant.VariantParams;
import com.rtg.variant.VariantParamsBuilder;
import com.rtg.variant.bayes.Description;
import com.rtg.variant.bayes.snp.ModelNoneFactory;
import com.rtg.variant.bayes.snp.ModelSnpFactory;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;

import junit.framework.TestCase;

/**
 */
public class IndividualSampleProcessorTest extends TestCase {

  private VariantAlignmentRecord getSamRecord(final int start) {
    final SAMRecord rec = new SAMRecord(new SAMFileHeader());
    rec.setReadString("TATT");
    rec.setAlignmentStart(start);
    rec.setCigarString("2M2D2M");
    return new VariantAlignmentRecord(rec);
  }

  public void test() throws IOException, InvalidParamsException {
    Diagnostic.setLogStream();
    final File output = FileUtils.createTempDir("testcheck", "variance_out");
    FileHelper.deleteAll(output);
    final File input = FileUtils.createTempDir("testcheck", "variance_in");
    try {
      final File samFile = new File(input, OUT_SAM);
      FileUtils.stringToFile(SAM1, samFile);
      final File templ = ReaderTestUtils.getDNADir(REF_SEQS);
      try {
        try (SequenceParams g = SequenceParams.builder().directory(templ).mode(SequenceMode.UNIDIRECTIONAL).create()) {
          final File outFile = new File("multisampleoutput");
          try {
            assertTrue(outFile.mkdir());
            final VariantParamsBuilder b = new VariantParamsBuilder();
            final ArrayList<File> list = new ArrayList<>();
            list.add(samFile);
            b.mapped(list);
            b.genome(g);
            b.outputParams(new OutputParams(outFile, false, false));
            final GenomePriorParams gpp = GenomePriorParams.builder().create();
            b.genomePriors(gpp);
            final VariantParams p = b.create();
            final ModelSnpFactory haploid = new ModelSnpFactory(p.genomePriors(), true);
            final ModelSnpFactory diploid = new ModelSnpFactory(p.genomePriors(), false);

            final byte[] template = new byte[10];
            final String name = "pufferfish";
            final SexMemo sexMemo = Utils.createSexMemo(p);
            final IndividualSampleFactory<Description> f = new IndividualSampleFactory<>(p, new DefaultMachineErrorChooser(), haploid, diploid, new ModelNoneFactory(), Sex.MALE, sexMemo);
            final IndividualSampleProcessor<?> task = f.make(name, template, 1, template.length);

            assertTrue(task.processAlignmentRecord(getSamRecord(1)));
            assertFalse(task.processAlignmentRecord(getSamRecord(-1)));
            assertFalse(task.processAlignmentRecord(getSamRecord(11)));
            assertTrue(task.processAlignmentRecord(getSamRecord(2)));
          } finally {
            assertTrue(FileHelper.deleteAll(outFile));
          }
        }
      } finally {
        assertTrue(FileHelper.deleteAll(templ));
      }
    } finally {
      assertTrue(FileHelper.deleteAll(output));
      assertTrue(FileHelper.deleteAll(input));
    }
  }
}
