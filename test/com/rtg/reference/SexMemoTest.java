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

package com.rtg.reference;

import static com.rtg.util.StringUtils.LS;
import static com.rtg.util.StringUtils.TAB;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;

import com.rtg.launcher.SequenceParams;
import com.rtg.util.test.RandomDna;
import com.rtg.mode.SequenceMode;
import com.rtg.reader.ReaderTestUtils;
import com.rtg.reader.SequencesReader;
import com.rtg.reference.ReferenceGenome.DefaultFallback;
import com.rtg.util.intervals.SequenceNameLocusSimple;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.io.TestDirectory;

import junit.framework.TestCase;

/**
 */
public class SexMemoTest extends TestCase {


  public void testDiploidDefault() throws IOException {
    Diagnostic.setLogStream();
    try (final TestDirectory genomeDir = new TestDirectory("sexmemo")) {
      ReaderTestUtils.getReaderDNA(">t\nacgt", genomeDir, null).close();
      final SequenceParams genomes = SequenceParams.builder().directory(genomeDir).mode(SequenceMode.UNIDIRECTIONAL).create();
      try (SequencesReader reader = genomes.reader()) {
        final SexMemo sx = new SexMemo(reader, DefaultFallback.DIPLOID);
        assertEquals(Ploidy.DIPLOID, sx.getEffectivePloidy(Sex.EITHER, "t"));
        assertEquals(Ploidy.NONE, sx.getEffectivePloidy(Sex.EITHER, "unknown"));
        assertEquals(Ploidy.DIPLOID, sx.getEffectivePloidy(Sex.MALE, "t"));
        assertEquals(Ploidy.DIPLOID, sx.getEffectivePloidy(Sex.FEMALE, "t"));
      }
    }
  }

  public void testHaploidDefault() throws IOException {
    Diagnostic.setLogStream();
    try (final TestDirectory genomeDir = new TestDirectory("sexmemo")) {
      ReaderTestUtils.getReaderDNA(">t\nacgt", genomeDir, null).close();
      final SequenceParams genomes = SequenceParams.builder().directory(genomeDir).mode(SequenceMode.UNIDIRECTIONAL).create();
      try (SequencesReader reader = genomes.reader()) {
        final SexMemo sx = new SexMemo(reader, DefaultFallback.HAPLOID);
        assertEquals(Ploidy.HAPLOID, sx.getEffectivePloidy(Sex.EITHER, "t"));
        assertEquals(Ploidy.NONE, sx.getEffectivePloidy(Sex.EITHER, "unknown"));
        assertEquals(Ploidy.HAPLOID, sx.getEffectivePloidy(Sex.MALE, "t"));
        assertEquals(Ploidy.HAPLOID, sx.getEffectivePloidy(Sex.FEMALE, "t"));
      }
    }
  }

  public void testRef() throws IOException {
    Diagnostic.setLogStream();
    try (final TestDirectory tempDir = new TestDirectory("sexmemo")) {
      final File genomeDir = ReaderTestUtils.getDNADir(">t\nacgt", tempDir, false, true, true);
      final SequenceParams genomes = SequenceParams.builder().directory(genomeDir).mode(SequenceMode.UNIDIRECTIONAL).create();
      try (SequencesReader reader = genomes.reader()) {
        final SexMemo sx = new SexMemo(reader, DefaultFallback.HAPLOID);
        assertEquals(Ploidy.DIPLOID, sx.getEffectivePloidy(Sex.EITHER, "t"));
        assertEquals(Ploidy.NONE, sx.getEffectivePloidy(Sex.EITHER, "unknown"));
        assertEquals(Ploidy.DIPLOID, sx.getEffectivePloidy(Sex.MALE, "t"));
        assertEquals(Ploidy.DIPLOID, sx.getEffectivePloidy(Sex.FEMALE, "t"));
      }
    }
  }

  public void testVariousPloids() throws IOException {
    //Diagnostic.setLogStream();
    try (final TestDirectory tempDir = new TestDirectory("sexmemo")) {
      final File genomeDir = ReaderTestUtils.getDNADir(">s1\n" + RandomDna.random(4000) + "\n>s2\n" + RandomDna.random(5000) + "\n>s3\nacgt\n", tempDir, false, true, true);
      final File refFile = new File(genomeDir, ReferenceGenome.REFERENCE_FILE);
      try (FileOutputStream fo = new FileOutputStream(refFile, true)) {
        fo.write(("female" + TAB + "seq" + TAB + "s1" + TAB + "diploid" + TAB + "circular" + LS
            + "female" + TAB + "seq" + TAB + "s2" + TAB + "none" + TAB + "circular" + LS
            + "male" + TAB + "seq" + TAB + "s1" + TAB + "haploid" + TAB + "circular" + TAB + "s2" + LS
            + "male" + TAB + "seq" + TAB + "s2" + TAB + "haploid" + TAB + "circular" + TAB + "s1" + LS
            + "male" + TAB + "dup" + TAB + "s1:1000-3000" + TAB + "s2:2000-4000" + LS
            + LS
            + "either" + TAB + "seq" + TAB + "s3" + TAB + "polyploid" + TAB + "circular" + LS
            + LS).getBytes());
      }
      final SequenceParams genomes = SequenceParams.builder().directory(genomeDir).mode(SequenceMode.UNIDIRECTIONAL).create();
      try (SequencesReader reader = genomes.reader()) {
        final SexMemo sx = new SexMemo(reader, DefaultFallback.HAPLOID);

        assertEquals(Ploidy.DIPLOID, sx.getRealPloidy(Sex.EITHER, "s1"));
        assertEquals(Ploidy.NONE, sx.getRealPloidy(Sex.EITHER, "unknown"));
        assertEquals(Ploidy.HAPLOID, sx.getRealPloidy(Sex.MALE, "s1"));
        assertEquals(Ploidy.DIPLOID, sx.getRealPloidy(Sex.FEMALE, "s1"));
        assertEquals(Ploidy.HAPLOID, sx.getRealPloidy(Sex.MALE, "s2"));
        assertEquals(Ploidy.NONE, sx.getRealPloidy(Sex.FEMALE, "s2"));
        assertEquals(Ploidy.HAPLOID, sx.getRealPloidy(Sex.MALE, "s2"));

        // Some PAR specific queries
        assertEquals(999, sx.getParBoundary(Sex.MALE, new SequenceNameLocusSimple("s1", 0, 2000)));
        assertEquals(3000, sx.getParBoundary(Sex.MALE, new SequenceNameLocusSimple("s1", 2000, 4000)));
        assertEquals(Ploidy.DIPLOID, sx.getRealPloidy(Sex.MALE, "s1", 1000));
        assertEquals(Ploidy.DIPLOID, sx.getEffectivePloidy(Sex.MALE, "s1", 1000));
        assertEquals(Ploidy.DIPLOID, sx.getRealPloidy(Sex.MALE, "s2", 2000));
        assertEquals(Ploidy.NONE, sx.getEffectivePloidy(Sex.MALE, "s2", 2000));

        assertEquals(Ploidy.POLYPLOID, sx.getRealPloidy(Sex.FEMALE, "s3"));
        assertEquals(Ploidy.HAPLOID, sx.getEffectivePloidy(Sex.FEMALE, "s3"));
      }
    }
  }

}
