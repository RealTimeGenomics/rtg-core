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
import com.rtg.reference.SexMemo;
import com.rtg.util.InvalidParamsException;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.io.FileUtils;
import com.rtg.util.test.FileHelper;
import com.rtg.variant.DefaultMachineErrorChooser;
import com.rtg.variant.GenomePriorParams;
import com.rtg.variant.VariantAlignmentRecord;
import com.rtg.variant.VariantParams;
import com.rtg.variant.VariantParamsBuilder;
import com.rtg.variant.bayes.Description;
import com.rtg.variant.bayes.NoAlleleBalance;
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
            b.genome(g.readerParams());
            b.outputParams(new OutputParams(outFile, false));
            final GenomePriorParams gpp = GenomePriorParams.builder().create();
            b.genomePriors(gpp);
            final VariantParams p = b.create();
            final ModelSnpFactory haploid = new ModelSnpFactory(p.genomePriors(), true, new NoAlleleBalance());
            final ModelSnpFactory diploid = new ModelSnpFactory(p.genomePriors(), false, new NoAlleleBalance());

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
