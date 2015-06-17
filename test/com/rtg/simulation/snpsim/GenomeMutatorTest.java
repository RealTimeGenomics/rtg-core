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
package com.rtg.simulation.snpsim;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;
import java.io.StringReader;
import java.util.Arrays;

import com.rtg.mode.DNA;
import com.rtg.mode.SequenceType;
import com.rtg.reader.PrereadType;
import com.rtg.reader.ReaderTestUtils;
import com.rtg.reader.Sdf2Fasta;
import com.rtg.reader.SdfWriter;
import com.rtg.reader.SequencesReader;
import com.rtg.reader.SequencesReaderFactory;
import com.rtg.reference.Ploidy;
import com.rtg.reference.Sex;
import com.rtg.simulation.snpsim.GenomeMutator.MutationComparator;
import com.rtg.simulation.snpsim.Mutation.DifferentMode;
import com.rtg.simulation.snpsim.Mutation.MutationType;
import com.rtg.util.InvalidParamsException;
import com.rtg.util.PortableRandom;
import com.rtg.util.StringUtils;
import com.rtg.util.TestUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.DiagnosticEvent;
import com.rtg.util.diagnostic.DiagnosticListener;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.IOUtils;
import com.rtg.util.io.TestDirectory;
import com.rtg.util.test.FileHelper;
import com.rtg.variant.GenomePriorParams;
import com.rtg.variant.VariantParams;
import com.rtg.variant.VariantParamsBuilder;
import com.rtg.vcf.VcfRecord;

import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;

/**
 * Tests corresponding class
 */
public class GenomeMutatorTest extends TestCase {

  private File mDir;
  private GenomePriorParams mPriors;
  static final String[] ERRSTR = new String[1];

  private static final VariantParams TEST_PARAMS = new VariantParamsBuilder().create();

  @Override
  public void setUp() throws IOException, InvalidParamsException {
    mDir = FileHelper.createTempDirectory();
    mPriors = GenomePriorParams.builder().genomePriors("testhumanprior").create();
  }

  @Override
  public void tearDown() {
    assertTrue(FileHelper.deleteAll(mDir));
    mPriors = null;
    mDir = null;
  }

  private static final String LS = StringUtils.LS;
  private static final String TB = "\t";
  private static final DiagnosticListener DIAG_LISTENER = new DiagnosticListener() {

    @Override
    public void handleDiagnosticEvent(final DiagnosticEvent<?> event) {
      ERRSTR[0] = event.getMessage();
    }

    @Override
    public void close() {
    }
  };

  private final class MyGenomeMutator extends GenomeMutator {

    MyGenomeMutator() throws InvalidParamsException, IOException {
      super(1, false, false, 1, TEST_PARAMS);
      final GenomePriorParams priors = GenomePriorParams.builder().genomePriors("human").create();
      setPriors(priors);
      mMutationGenerator = new GenomeMutatorPriors(mPriors, Ploidy.DIPLOID);

      final Mutation mOmo = new Mutation(22, Mutation.MutationType.SNP, false,
          DifferentMode.ONE_ONLY, GenDiffMode.TWIN_ONLY, 1, 1, null, null);
      final Mutation mEtero = new Mutation(22, MutationType.SNP, false,
          DifferentMode.ONE_ONLY, GenDiffMode.TWIN_ONLY, 1, 1, null, null);
      final String[] strEtero = {"A", "C", "G"};
      final String[] strOmo = {"A", "C", "G"};
      try {
        mutationLine(null, null, null, mOmo, 'N', false, strOmo);
        fail();
      } catch (final IllegalArgumentException e) {
        assertEquals("Unpossible", e.getMessage());
      }
      try {
        mutationLine(null, null, "a", null, 'N', false, strOmo);
        fail();
      } catch (final IllegalArgumentException e) {
        assertEquals("Unpossible", e.getMessage());
      }
      try {
        mutationLine(null, null, "a", mOmo, 'N', false, "A");
        fail();
      } catch (final IllegalArgumentException e) {
        assertEquals("Unpossible", e.getMessage());
      }
      try {
        mutationLine(null, null, "a", mOmo, 'N', false, "A", "A", "A", "A");
        fail();
      } catch (final IllegalArgumentException e) {
        assertEquals("Unpossible", e.getMessage());
      }
      try {
        mutationLine(null, null, "a", mOmo, 'N', false, strEtero);
        fail();
      } catch (final IllegalArgumentException e) {
        assertEquals("Unpossible", e.getMessage());
      }
      try {
        mutationLine(null, null, "a", mEtero, 'N', false, strOmo);
        fail();
      } catch (final IllegalArgumentException e) {
        assertEquals("Unpossible", e.getMessage());
      }
      try {
        dnasToString(null);
      } catch (final IllegalArgumentException e) {
        assertEquals("MNP dnas shouldnt be null", e.getMessage());
      }
      for (int i = 0; i < 10; i++) {
        final DNA[] dna = new DNA[1];
        dna[0] = randomDNA(new PortableRandom(2));
        assertTrue(dnasToString(dna).matches("[ACGT]"));
      }
    }
  }

  private static byte[] dnaToSeq(String dnaString) {
    final byte[] seq = new byte[dnaString.length()];
    for (int i = 0; i < seq.length; i++) {
      seq[i] = (byte) DNA.valueOf(Character.toString(dnaString.charAt(i))).ordinal();
    }
    return seq;
  }

  private static class MyGenomeMutator2 extends GenomeMutator {
    protected MyGenomeMutator2(GenomePriorParams priors) throws InvalidParamsException {
      super(22, false, true, 1, TEST_PARAMS);
      mMutationGenerator = new GenomeMutatorPriors(priors, Ploidy.DIPLOID);
      final byte[] seq = dnaToSeq("ACCCTTGGG");
      DNA[] dna = getMutatedMnp(seq, getDNAFromBytes(seq), 0, seq.length, seq.length);
      assertEquals("CCCTTGG", dnasToString(dna));
      dna = getMutatedMnp(seq, getDNAFromBytes(seq), 0, seq.length, seq.length);
      assertEquals("CCCCTTGGA", dnasToString(dna));
      dna = getMutatedMnp(seq, getDNAFromBytes(seq), 0, seq.length, seq.length);
      assertEquals("CCCTTGG", dnasToString(dna));
      moreTests();

    }

    protected MyGenomeMutator2(final int seed, final boolean verbose, boolean simpleMnps) {
      super(seed, verbose, simpleMnps, 1, TEST_PARAMS);
    }

    private void moreTests() {
      //snp and mnp
      final Mutation m0 = new Mutation(4, MutationType.SNP, false, DifferentMode.HOMOZYGOUS,
          GenDiffMode.BOTH_SAME, 1, 1, null, null);
      assertTrue(MutationTest.integrity(m0));
      check(m0, 4, 1, 0, "AAACxCTTT", true, "seq1" + TB + "5" + TB + "." + TB + "C" + TB + "A" + TB + "." + TB + "PASS" + TB + "." + TB + "GT:GQ" + TB + "1/1:3");
      check(m0, 4, 1, 0, "AAACxCTTT", false, "seq1" + TB + "5" + TB + "." + TB + "C" + TB + "G" + TB + "." + TB + "PASS" + TB + "." + TB + "GT:GQ" + TB + "1/1:3");
      final Mutation m1 = new Mutation(4, MutationType.MNP, false, DifferentMode.HOMOZYGOUS,
          GenDiffMode.BOTH_SAME, 3, 3, null, null);
      assertTrue(MutationTest.integrity(m1));
      checkMnp(m1, 4, "AAACxxTT", true);
      checkMnp(m1, 4, "AAACxxxTT", false);
      final Mutation m3 = new Mutation(6, MutationType.INSERT, false, DifferentMode.HOMOZYGOUS,
          GenDiffMode.BOTH_SAME, 2, 2, null, null);
      assertTrue(MutationTest.integrity(m3));
      check(m3, -1, -1, 2, "AAACCCxxTTT", true, "seq1" + TB + "6" + TB + "." + TB + "C" + TB + "CGT" + TB + "." + TB + "PASS" + TB + "." + TB + "GT:GQ" + TB + "1/1:3");
      check(m3, -1, -1, 2, "AAACCCxxTTT", false, "seq1" + TB + "6" + TB + "." + TB + "C" + TB + "CCG" + TB + "." + TB + "PASS" + TB + "." + TB + "GT:GQ" + TB + "1/1:3");
      final Mutation m4 = new Mutation(6, MutationType.DELETE, false, DifferentMode.HOMOZYGOUS,
          GenDiffMode.BOTH_SAME, 2, 2, null, null);
      check(m4, 7, 2, -2, "AAACCCT", true, "seq1" + TB + "6" + TB + "." + TB + "CTT" + TB + "C" + TB + "." + TB + "PASS" + TB + "." + TB + "GT:GQ" + TB + "1/1:3");
      check(m4, 7, 2, -2, "AAACCCT", false, "seq1" + TB + "6" + TB + "." + TB + "CTT" + TB + "C" + TB + "." + TB + "PASS" + TB + "." + TB + "GT:GQ" + TB + "1/1:3");
      assertTrue(MutationTest.integrity(m4));

      final Mutation m5 = new Mutation(4, MutationType.SNP, true, DifferentMode.ONE_ONLY,
          GenDiffMode.FIRST_ONLY, 1, 0, null, null);
      assertTrue(MutationTest.integrity(m5));
      check(m5, 4, 1, 0, "AAACxCTTT", true, "seq1" + TB + "5" + TB + "." + TB + "C" + TB + "T" + TB + "." + TB + "PASS" + TB + "." + TB + "GT:GQ" + TB + "1/0:3");
      check(m5, -1, 0, 0, "AAACCCTTT", false, "seq1" + TB + "5" + TB + "." + TB + "C" + TB + "G" + TB + "." + TB + "PASS" + TB + "." + TB + "GT:GQ" + TB + "1/0:3");
      final Mutation m5a = new Mutation(4, MutationType.SNP, true, DifferentMode.ONE_ONLY,
          GenDiffMode.TWIN_ONLY, 0, 1, null, null);
      //
      assertTrue(MutationTest.integrity(m5a));
      check(m5a, -1, 0, 0, "AAACCCTTT", true, "seq1" + TB + "5" + TB + "." + TB + "C" + TB + "T" + TB + "." + TB + "PASS" + TB + "." + TB + "GT:GQ" + TB + "0/1:3");
      check(m5a, 4, 1, 0, "AAACxCTTT", false, "seq1" + TB + "5" + TB + "." + TB + "C" + TB + "G" + TB + "." + TB + "PASS" + TB + "." + TB + "GT:GQ" + TB + "0/1:3");
      final Mutation m6 = new Mutation(4, MutationType.MNP, true, DifferentMode.DIFFERENT,
          GenDiffMode.DIFFERENT, 3, 2, null, null);
      assertTrue(MutationTest.integrity(m6));
      checkMnp(m6, 4, "AAACxxxTT", true);
      checkMnp(m6, 4, "AAACxxTTT", false);

      // insert
      Mutation m7 = new Mutation(5, MutationType.INSERT, true, DifferentMode.DIFFERENT,
          GenDiffMode.DIFFERENT, 1, 2, null, null);
      assertTrue(MutationTest.integrity(m7));
      check(m7, -6, -1, 1, "AAACCxCTTT", true, "seq1" + TB + "5" + TB + "." + TB + "C" + TB + "CG,CGC" + TB + "." + TB + "PASS" + TB + "." + TB + "GT:GQ" + TB + "1/2:3");

      check(m7, -6, -1, 2, "AAACCxxCTTT", false, "seq1" + TB + "5" + TB + "." + TB + "C" + TB + "CA,CGT" + TB + "." + TB + "PASS" + TB + "." + TB + "GT:GQ" + TB + "1/2:3");

      m7 = new Mutation(7, MutationType.INSERT, false, DifferentMode.HOMOZYGOUS,
          GenDiffMode.BOTH_SAME, 2, 2, null, null);
      assertTrue(MutationTest.integrity(m7));
      check(m7, -8, -2, 2, "AAACCCTxxTT", true, "seq1" + TB + "7" + TB + "." + TB + "T" + TB + "TGT" + TB + "." + TB + "PASS" + TB + "." + TB + "GT:GQ" + TB + "1/1:3");

      check(m7, -8, -2, 2, "AAACCCTxxTT", false, "seq1" + TB + "7" + TB + "." + TB + "T" + TB + "TCT" + TB + "." + TB + "PASS" + TB + "." + TB + "GT:GQ" + TB + "1/1:3");

      m7 = new Mutation(5, MutationType.INSERT, true, DifferentMode.ONE_ONLY,
          GenDiffMode.TWIN_ONLY, 0, 2, null, null);
      assertTrue(MutationTest.integrity(m7));
      check(m7, -6, 0, 0, "AAACCCTTT", true, "seq1" + TB + "5" + TB + "." + TB + "C" + TB + "CTA" + TB + "." + TB + "PASS" + TB + "." + TB + "GT:GQ" + TB + "0/1:3");

      check(m7, -6, -1, 2, "AAACCxxCTTT", false, "seq1" + TB + "5" + TB + "." + TB + "C" + TB + "CGC" + TB + "." + TB + "PASS" + TB + "." + TB + "GT:GQ" + TB + "0/1:3");

      m7 = new Mutation(5, MutationType.INSERT, true, DifferentMode.ONE_ONLY,
          GenDiffMode.FIRST_ONLY, 1, 0, null, null);
      assertTrue(MutationTest.integrity(m7));
      check(m7, -6, -1, 1, "AAACCxCTTT", true, "seq1" + TB + "5" + TB + "." + TB + "C" + TB + "CT" + TB + "." + TB + "PASS" + TB + "." + TB + "GT:GQ" + TB + "1/0:3");
      check(m7, -6, -1, 0, "AAACCCTTT", false, "seq1" + TB + "5" + TB + "." + TB + "C" + TB + "CG" + TB + "." + TB + "PASS" + TB + "." + TB + "GT:GQ" + TB + "1/0:3");

      // delete
      final byte[] genome = getBytesFromDNA(GENOMESTR);

      Mutation m8 = new Mutation(1, MutationType.DELETE, true, DifferentMode.DIFFERENT,
          GenDiffMode.DIFFERENT, 3, 6, null, null);
      final Mutation m9 = new Mutation(8, MutationType.SNP, false, DifferentMode.HOMOZYGOUS,
          GenDiffMode.BOTH_SAME, 1, 1, null, null);

      assertTrue(MutationTest.integrity(m8));
      check(m8, -1, -1, -3, "AxxxTT", true, "seq1" + TB + "1" + TB + "." + TB + "AAACCCT" + TB + "ACCT,A" + TB + "." + TB + "PASS" + TB + "." + TB + "GT:GQ" + TB + "1/2:3");
      check(m8, 1, 8, -6, "ATT", false, "seq1" + TB + "1" + TB + "." + TB + "AAACCCT" + TB + "ACCT,A" + TB + "." + TB + "PASS" + TB + "." + TB + "GT:GQ" + TB + "1/2:3");
      assertTrue(validPos(new Mutation[] {m9}, 1, m8, genome));
      m8 = new Mutation(1, MutationType.DELETE, false, DifferentMode.HOMOZYGOUS,
          GenDiffMode.BOTH_SAME, 3, 3, null, null);
      assertTrue(MutationTest.integrity(m8));
      check(m8, -1, -1, -3, "AxxxTT", true, "seq1" + TB + "1" + TB + "." + TB + "AAAC" + TB + "A" + TB + "." + TB + "PASS" + TB + "." + TB + "GT:GQ" + TB + "1/1:3");
      check(m8, 1, 8, -3, "AxxxTT", false, "seq1" + TB + "1" + TB + "." + TB + "AAAC" + TB + "A" + TB + "." + TB + "PASS" + TB + "." + TB + "GT:GQ" + TB + "1/1:3");
      assertTrue(validPos(new Mutation[] {m9}, 1, m8, genome));
      m8 = new Mutation(1, MutationType.DELETE, true, DifferentMode.ONE_ONLY,
          GenDiffMode.TWIN_ONLY, 0, 6, null, null);
      assertTrue(MutationTest.integrity(m8));
      check(m8, -1, -1, 0, "AxxxxxxTT", true, "seq1" + TB + "1" + TB + "." + TB + "AAACCCT" + TB + "A" + TB + "." + TB + "PASS" + TB + "." + TB + "GT:GQ" + TB + "0/1:3");
      check(m8, 1, 8, -6, "ATT", false, "seq1" + TB + "1" + TB + "." + TB + "AAACCCT" + TB + "A" + TB + "." + TB + "PASS" + TB + "." + TB + "GT:GQ" + TB + "0/1:3");
      assertTrue(validPos(new Mutation[] {m9}, 1, m8, genome));
      m8 = new Mutation(1, MutationType.DELETE, true, DifferentMode.ONE_ONLY,
          GenDiffMode.FIRST_ONLY, 3, 0, null, null);
      assertTrue(MutationTest.integrity(m8));
      check(m8, -1, -1, -3, "AxxxTT", true, "seq1" + TB + "1" + TB + "." + TB + "AAAC" + TB + "A" + TB + "." + TB + "PASS" + TB + "." + TB + "GT:GQ" + TB + "1/0:3");
      check(m8, 1, -1, 0, "AxxxCCTT", false, "seq1" + TB + "1" + TB + "." + TB + "AAAC" + TB + "A" + TB + "." + TB + "PASS" + TB + "." + TB + "GT:GQ" + TB + "1/0:3");
      assertTrue(validPos(new Mutation[] {m9}, 1, m8, genome));

      final byte[] nsGenome = getBytesFromDNA("NNNNNNNN");
      assertFalse(validPos(new Mutation[] {m0}, 0, m8, nsGenome));


      // test first mnp shorter than second
      final Mutation m10 = new Mutation(4, MutationType.MNP, true, DifferentMode.DIFFERENT,
          GenDiffMode.DIFFERENT, 2, 3, null, null);
      assertTrue(MutationTest.integrity(m6));
      checkMnp(m10, 4, "AAACxxTTT", true);
      checkMnp(m10, 4, "AAACxxTTT", false);

      final Mutation m11 = new Mutation(4, MutationType.MNP, false, DifferentMode.HOMOZYGOUS,
          GenDiffMode.BOTH_SAME, 3, 3, null, null);
      assertTrue(MutationTest.integrity(m11));
      checkMnp(m11, 4, "AAACxxxTT", true);
      checkMnp(m11, 4, "AAACxxxTT", false);

    }
    static final String GENOMESTR = "AAACCCTTT";

    private String mMutantStr = "";

    @Override
    VcfRecord mutationLine(SdfWriter h1, SdfWriter h2, String seqName, Mutation m, char prevRefNt, boolean verbose, String... newAndOld) throws IOException {
      VcfRecord rec = super.mutationLine(h1, h2, seqName, m, prevRefNt, verbose, newAndOld);
      mMutantStr = rec.toString();
      return rec;
    }

    private void check(Mutation m1, final int pos, final int numMismatch, final int increase,
        String expectedGenome, boolean testFirst, String mapFileExpected) {
      try {
        try (final TestDirectory temp = new TestDirectory()) {
          final File mutant = new File(temp, "mutant2");
          final File twinMutant = new File(temp, "twinmutant");
          final byte[] genome = getBytesFromDNA(GENOMESTR);
          mMutantStr = "";

          try (SdfWriter output = new SdfWriter(mutant, 100, PrereadType.UNKNOWN, false, true, false, SequenceType.DNA); SdfWriter twinOutput = new SdfWriter(twinMutant, 100, PrereadType.UNKNOWN, false, true, false, SequenceType.DNA)) {
            output.startSequence("seq1");
            twinOutput.startSequence("seq1");

            mutateSequence(output, twinOutput, genome, new Mutation[]{m1}, "seq1");
            output.endSequence();
            twinOutput.endSequence();
          }
          try (SequencesReader t = SequencesReaderFactory.createDefaultSequencesReader(twinMutant)) {
            t.read(0);
            //System.err.println("Twi: " + Arrays.toString(getDNAFromBytes(twin)));
          }
          try (SequencesReader r = testFirst ? SequencesReaderFactory.createDefaultSequencesReader(mutant) : SequencesReaderFactory.createDefaultSequencesReader(twinMutant)) {
            final byte[] newMutant = r.read(0);
            int mutantPos;
            if (pos >= 0) {
              mutantPos = findMismatchPosition(genome, newMutant);
              if (mutantPos != -1) {
                assertEquals(pos, mutantPos);
              }
            }
            assertEquals(genome.length + increase, newMutant.length);
            if (numMismatch >= 0) {
              assertEquals(numMismatch, numMismatch(genome, newMutant));
            }
            if (expectedGenome.length() > 0) {
              assertTrue(checkEqualGenome(expectedGenome, getDNAStringFromBytes(newMutant)));
            }
            TestUtils.containsAll(mMutantStr, mapFileExpected);
          }
        }
      } catch (final IOException e) {
        fail(e.getMessage());
      }
    }

    protected void checkMnp(Mutation m1Mnp, final int posMnp,
                            String expectedGenome, boolean testFirst) {
      try {
        try (final TestDirectory temp = new TestDirectory()) {
          final File mutant = new File(temp, "mutant2");
          final File twinMutant = new File(temp, "twinmutant");
          final String genomeStr = "AAACCCTTT";
          final byte[] genome = getBytesFromDNA(genomeStr);
          mMutantStr = "";

          try (SdfWriter twinOutput = new SdfWriter(twinMutant, 100, PrereadType.UNKNOWN, false, true, false, SequenceType.DNA); SdfWriter output = new SdfWriter(mutant, 100, PrereadType.UNKNOWN, false, true, false, SequenceType.DNA)) {
            output.startSequence("seq1");
            twinOutput.startSequence("seq1");

            mutateSequence(output, twinOutput, genome, new Mutation[]{m1Mnp}, "seq1");
            output.endSequence();
            twinOutput.endSequence();
          }
          //System.err.println(FileUtils.fileToString(mapping));

          try (SequencesReader r = testFirst ? SequencesReaderFactory.createDefaultSequencesReader(mutant) : SequencesReaderFactory.createDefaultSequencesReader(twinMutant)) {
            final byte[] newMutant = r.read(0);
            //            System.err.println(Arrays.toString(genome));
            //            System.err.println(Arrays.toString(b));
            //System.err.println("Exp: " + expectedGenome + "\nWas: " + getDNAFromBytes(newMutant));
            if (posMnp >= 0) {
              final int mismatch = findMismatchPosition(genome, newMutant);
              assertTrue(mismatch >= posMnp || mismatch == -1);
            }
            if (expectedGenome.length() > 0) {
              assertTrue(checkEqualGenomeStart(expectedGenome, getDNAStringFromBytes(newMutant)));
              assertTrue(checkEqualGenomeEnd(expectedGenome, getDNAStringFromBytes(newMutant)));
            }
          }
        }
      } catch (final IOException e) {
        //System.err.println(e.getMessage());
        fail(e.getMessage());
      }
    }

    private int findMismatchPosition(byte[] genome, byte[] mutant) {
      for (int i = 0; i < genome.length; i++) {
        if (mutant.length <= i) {
          return i;
        }
        if (genome[i] != mutant[i]) {
          return i;
        }
      }
      return -1;
    }

    private boolean checkEqualGenome(String genome, String mutant) {
      for (int i = 0; i < genome.length(); i++) {
        final String c = genome.substring(i, i + 1);
        if (c.compareTo("x") != 0) {
          if (c.compareTo(mutant.substring(i, i + 1)) != 0) {
            return false;
          }
        }
      }
      return true;
    }

    private boolean checkEqualGenomeStart(String genome, String mutant) {
      for (int i = 0; i < genome.length(); i++) {
        final String c = genome.substring(i, i + 1);
        if (c.compareTo("x") != 0) {
          if (c.compareTo(mutant.substring(i, i + 1)) != 0) {
            return false;
          }
        } else {
          return true;
        }
      }
      return true;

    }

    private boolean checkEqualGenomeEnd(String genome, String mutant) {
      for (int i = 0; i > genome.length(); i++) {
        final String c = genome.substring(genome.length() - 1 - i, genome.length() - i);
        if (c.compareTo("x") != 0) {
          if (c.compareTo(mutant.substring(mutant.length() - 1 - i, mutant.length() - i)) != 0) {
            return false;
          }
        } else {
          return true;
        }
      }
      return true;

    }

    private int numMismatch(byte[] genome, byte[] mutant) {
      int num = 0;
      for (int i = 0; i < genome.length; i++) {
        if (mutant.length <= i) {
          num++;
        } else if (genome[i] != mutant[i]) {
          num++;
        }
      }
      return num;
    }

    private byte[] getBytesFromDNA(String dna) {
      final byte[] b = new byte[dna.length()];
      for (int i = 0; i < dna.length(); i++) {
        b[i] = (byte) DNA.valueOf("" + dna.charAt(i)).ordinal();
      }
      return b;
    }

    private String getDNAStringFromBytes(byte[] b) {
      final StringBuilder str = new StringBuilder();
      for (final byte v : b) {
        str.append(DNA.values()[v].toString());
      }
      return str.toString();
    }

    private DNA[] getDNAFromBytes(byte[] b) {
      final DNA[] dna = new DNA[b.length];
      for (int i = 0; i < b.length; i++) {
        dna[i] = DNA.values()[b[i]];
      }
      return dna;
    }
  }

  private final class MyGenomeMutator3 extends MyGenomeMutator2 {
    MyGenomeMutator3() throws InvalidParamsException {
      super(22, false, false);
      mMutationGenerator = new GenomeMutatorPriors(mPriors, Ploidy.DIPLOID);
      // test first mnp shorter than second
      final Mutation m10 = new Mutation(4, MutationType.MNP, true, DifferentMode.DIFFERENT,
          GenDiffMode.DIFFERENT, 2, 3, null, null);
      assertTrue(MutationTest.integrity(m10));
      checkMnp(m10, 4, "AAACxxTTT", true);
      checkMnp(m10, 4, "AAACxxTTT", false);
    }
  }

  public void testMyGenomeMutator() throws IOException, InvalidParamsException {
    Diagnostic.setLogStream();
    new MyGenomeMutator();
  }

  public void testMyGenomeMutator2() throws IOException, InvalidParamsException {
    Diagnostic.setLogStream();
    new MyGenomeMutator2(mPriors);
  }

  public void testMyGenomeMutator3() throws IOException, InvalidParamsException {
    Diagnostic.setLogStream();
    new MyGenomeMutator3();
  }

  private static class MyMutationComparator extends MutationComparator {

    MyMutationComparator() {
      final Mutation m1 = new Mutation(1, MutationType.SNP, true,
          DifferentMode.ONE_ONLY, GenDiffMode.FIRST_ONLY, 2, 0,
        null, null);
      final Mutation m2 = new Mutation(2, MutationType.SNP, true,
          DifferentMode.ONE_ONLY, GenDiffMode.FIRST_ONLY, 2, 0,
        null, null);
      assertTrue(compare(m1, m2) < 0);
    }
  }

  public static void main(final String[] args) {
    junit.textui.TestRunner.run(suite());
  }

  public static Test suite() {
    final TestSuite suite = new TestSuite();

    suite.addTestSuite(GenomeMutatorTest.class);
    return suite;
  }
  //private static final OutputStream NULL_STREAM = SimpleTestUtils.getNullOutputStream();
  private static final PrintStream NULL_PRINTSTREAM = TestUtils.getNullPrintStream();

  public void testSnpMutator() throws IOException {
    final File temp = FileHelper.createTempDirectory();
    try {
      final String genomeTides = "ACGTACGATCAGCATCTGACATGCTAACGGTCATCGCGGCATTTACGGCACGTCATAGTCGCATGCGGCATTTATAGCGGCGCGTAAATCGTCATGATTGCAGCTAGCATCGTGTCACA";

      final String genome = ">a\n" + genomeTides + "\n";
      final File in = ReaderTestUtils.getDNASubDir(genome, mDir);
      try {
        try (SequencesReader dsr = SequencesReaderFactory.createDefaultSequencesReader(in)) {
          final GenomeMutator gm = new GenomeMutator(1, false, false, 1, TEST_PARAMS);
          final File mutant = new File(temp, "mutant");
          final File mutations = new File(temp, "mutations.mapping");
          gm.setRates(1.0, 0.0, 0.0);
          try (OutputStream os = new FileOutputStream(mutations)) {
            gm.mutateCount(dsr, Sex.EITHER, mutant, null, os, 40);
          }
          toFasta(mutant);
          final File output = new File(mutant.getPath() + ".fasta");
          final String outStr = IOUtils.readAll(output);
          //System.err.println(outStr);
          //System.err.println("123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890");
          final String stripped = outStr.replaceAll(">[^\n]*\n", "").replaceAll("\\s", "");
          final String mapStr = IOUtils.readAll(mutations);
          //assertTrue(!genomeTides.equals(stripped));
          final char[] mod = stripped.toCharArray();
          final BufferedReader sr = new BufferedReader(new StringReader(mapStr));
          String line;
          while ((line = sr.readLine()) != null) {
            if (!line.substring(0, 1).equals("#")) {
              final String[] parts = line.split("\t");
              assertEquals("a", parts[0]);
              final int pos = Integer.parseInt(parts[1]) - 1;
              //System.err.println("" + pos + ": " + line);

              assertEquals(".", parts[2]);
              assertEquals(parts[3].charAt(0), genomeTides.charAt(pos));
              assertEquals(parts[4].charAt(0), stripped.charAt(pos));
              mod[pos] = parts[3].charAt(0);
              assertEquals("1:3", parts[9]);
            }
          }
          assertTrue(stripped.matches("^[ACGT]+$"));
          assertEquals(genomeTides, new String(mod));
          //assertTrue(count > 0);
          TestUtils.containsAll(gm.statisticsLines(),
            "Passed Filters               : 40",
            "Failed Filters               : 0",
            "SNPs                         : 40",
            "Haploid SNPs                 : 40"
          );
        }
      } finally {
        assertTrue(FileHelper.deleteAll(in));
      }
    } finally {
      assertTrue(FileHelper.deleteAll(temp));
    }
  }

  public void testTooHighCount() throws IOException {
    final String[] errStr = new String[1];

    final DiagnosticListener dl = new DiagnosticListener() {

      @Override
      public void handleDiagnosticEvent(final DiagnosticEvent<?> event) {
        errStr[0] = event.getMessage();
      }

      @Override
      public void close() {
      }
    };
    Diagnostic.addListener(dl);
    try {
      Diagnostic.setLogStream();
      final File temp = FileUtils.createTempDir("genomemutatortest", "toohighcount");
      final String genome = ""
          + ">a" + LS
          + "ACGTACG" + LS;
      final File in = ReaderTestUtils.getDNASubDir(genome, mDir);
      try {
        final SequencesReader dsr = SequencesReaderFactory.createDefaultSequencesReader(in);
        final GenomeMutator gm = new GenomeMutator(1, false, false, 1, TEST_PARAMS);
        final File mutations = new File(temp, "mutations.mapping");
        final File mutant = new File(temp, "mutant");
        try (OutputStream s = new FileOutputStream(mutations)) {
          assertEquals(1, gm.mutateCount(dsr, Sex.EITHER, mutant, null, s, 20));
        }
        TestUtils.containsAll(errStr[0], "Count gives a mutation rate of greater than 1 mutation per nucleotide");
      } finally {
        assertTrue(FileHelper.deleteAll(temp));
      }
    } finally {
      Diagnostic.removeListener(dl);
    }
  }

  public void testTooSmallRate() throws IOException {
    //final String[] errStr = new String[1];

    Diagnostic.addListener(DIAG_LISTENER);
    try {
      Diagnostic.setLogStream();
      final String genome = ""
          + ">a" + LS
          + "ACGTGGG" + LS;
      final File in = ReaderTestUtils.getDNASubDir(genome, mDir);
      try (SequencesReader dsr = SequencesReaderFactory.createDefaultSequencesReader(in)) {
        final GenomeMutator gm = new GenomeMutator(1, false, false, 1, TEST_PARAMS);
        final File mutations = new File(mDir, "mutations.mapping");
        final File mutant = new File(mDir, "mutant");
        try (OutputStream s = new FileOutputStream(mutations)) {
          assertEquals(0, gm.mutateRate(dsr, Sex.EITHER, mutant, null, s, 0.1));
        }
        TestUtils.containsAll(ERRSTR[0], "Mutation rate is less than 1 mutation per sequence");
      }
    } finally {
      Diagnostic.removeListener(DIAG_LISTENER);
    }
  }

  public void testMultipleSequences() throws IOException {
    final File temp = FileHelper.createTempDirectory();
    try {
      final String genome = ""
          + ">a\n"
          + "ACGTACGATCAGCATCTGACATGCTAACGGTCATCGCGGCATTTACGGCACGTCATAGTCGCATGCGGCATTTATAGCGGCGCGTAAATCGT" + "\n"
          + ">b\n"
          + "CATGATTGCAGCTAGCATCGTGTCACA" + "\n";
      final File in = ReaderTestUtils.getDNASubDir(genome, mDir);
      try {
        final SequencesReader dsr = SequencesReaderFactory.createDefaultSequencesReader(in);
        final GenomeMutator gm = new GenomeMutator(1, false, false, 1, TEST_PARAMS);
        gm.setRates(1.0, 0.0, 0.0);
        final File mutations = new File(temp, "mutations.mapping");
        final File mutant = new File(temp, "mutant");
        try (OutputStream os = new FileOutputStream(mutations)) {
          gm.mutateCount(dsr, Sex.EITHER, mutant, null, os, 5);
        }
        toFasta(mutant);
        final BufferedReader r1 = new BufferedReader(new FileReader(mutant.getPath() + ".fasta"));
        final BufferedReader r2 = new BufferedReader(new StringReader(genome));
        final int count = compareCount(r1, r2);
        final BufferedReader r = new BufferedReader(new FileReader(mutations));
        final int lcnt = lineCount(r);
        assertEquals(count, lcnt);
        r1.close();
        r2.close();
        r.close();
        dsr.close();
      } finally {
        assertTrue(FileHelper.deleteAll(in));
      }
    } finally {
      assertTrue(FileHelper.deleteAll(temp));
    }
  }

  private void toFasta(final File preread) throws IOException {
    new Sdf2Fasta().mainInit(new String[] {"-i", preread.toString(), "-Z", "-o", preread.getPath() + ".fasta"}, FileUtils.getStdoutAsOutputStream(), System.out);
  }

  private int compareCount(final BufferedReader r1, final BufferedReader r2) throws IOException {
    int count = 0;
    while (true) {
      final String line1 = r1.readLine();
      final String line2 = r2.readLine();
      if (line1 == null && line2 == null) {
        break;
      }
      assertNotNull(line1);
      assertNotNull(line2);
      assertEquals(line1.length(), line2.length());
      for (int i = 0; i < line1.length(); i++) {
        if (line1.charAt(i) != line2.charAt(i)) {
          count++;
        }
      }
    }
    return count;
  }

  private int lineCount(File f) throws IOException {
    try (BufferedReader br = new BufferedReader(new FileReader(f))) {
      return lineCount(br);
    }
  }
  private int lineCount(final BufferedReader r) throws IOException {
    int count = 0;
    String line;
    while ((line = r.readLine()) != null) {
      if (!line.substring(0, 1).equals("#")) {
        //System.err.println(line);
        count++;
      }
    }
    return count;
  }

  public void testMnpMutationDiploid() throws IOException {
    final File temp = FileHelper.createTempDirectory();
    try {
      final String genomeTides = "ACGTACGATCAGCATCTGACATCGGCTACTACGGCATTGCAATCGGCTACGAT";
      final String genome = ">a\n" + genomeTides + "\n";
      final File in = ReaderTestUtils.getDNASubDir(genome, mDir);
      try {
        try (SequencesReader dsr = SequencesReaderFactory.createDefaultSequencesReader(in)) {
          final File mutations = new File(temp, "mutations.mapping");
          final File mutant = new File(temp, "mutant");
          final GenomeMutator gm = new GenomeMutator(1, false, false, 1, TEST_PARAMS);
          gm.setRates(0.0, 1.0, 0.0);
          final OutputStream os = new FileOutputStream(mutations);
          try {
            gm.mutateCount(dsr, Sex.EITHER, mutant, null, os, 2);
          } finally {
            os.close();
          }
          toFasta(mutant);
          final File output = new File(mutant.getPath() + ".fasta");
          final String outStr = IOUtils.readAll(output);
          final String stripped = outStr.replaceAll(">[^\n]*\n", "").replaceAll("\\s", "");
          final String mapStr = IOUtils.readAll(mutations);
          //assertTrue(!genomeTides.equals(stripped));
          final BufferedReader sr = new BufferedReader(new StringReader(mapStr));
          String line;
          int posDiff = 0;
          while ((line = sr.readLine()) != null) {
            if (!line.substring(0, 1).equals("#")) {
              final String[] parts = line.split("\t");
              assertEquals("a", parts[0]);
              final int pos = Integer.parseInt(parts[1]) - 1;
              assertEquals(".", parts[2]);
              assertEquals("1:3", parts[9]);
              assertEquals(parts[3], genomeTides.substring(pos, pos + parts[3].length()));
              assertEquals(parts[4], stripped.substring(pos + posDiff, pos + posDiff + parts[4].length()));
              posDiff += parts[4].length() - parts[3].length();
              //assertEquals(parts[3].length(), parts[4].length());
            }
          }
          assertTrue(stripped.matches("^[ACGT]+$"));
          //mutator will skip mutations if it has already moved past mutation point
          //assertTrue(count > 0);

        }
      } finally {
        assertTrue(FileHelper.deleteAll(in));
      }
    } finally {
      assertTrue(FileHelper.deleteAll(temp));
    }
  }

  public void testIndelMutation() throws IOException {
    final File temp = FileHelper.createTempDirectory();
    try {
      //final String genomeTides = "ACGTACGATCAGCATCTGACATCGGCTACTACGGCATTGCAATCGGCTACGATATGAGTCTACGCTAATGCAGAGCTACGTGTCTAGACTCAGTCAGGTCGCATGACTCATATGAGTAGTACGGCATGCACGATC";
      final String genomeTides = "TCTACGCTTAGTACGGCATGCACGATC";
      //final String xxnums      = "012345678901234567890123456";
      final String genome = ">sequence1\n" + genomeTides + "\n";
      final File in = ReaderTestUtils.getDNASubDir(genome, mDir);
      try {
        for (int s = 1; s < 5; s++) {
          try (SequencesReader dsr = SequencesReaderFactory.createDefaultSequencesReader(in)) {
            //System.err.println("seed " + s);
            final File mutations = new File(temp, "mutations.mapping");
            final File mutant = new File(temp, "mutant" + s);
            final GenomeMutator gm = new GenomeMutator(s, false, false, 1, TEST_PARAMS);
            gm.setRates(0.0, 0.0, 1.0);

            gm.mutateCount(dsr, Sex.EITHER, mutant, null, new FileOutputStream(mutations), 1);
            toFasta(mutant);
            final File output = new File(mutant.getPath() + ".fasta");
            final String outStr = IOUtils.readAll(output);
            final String stripped = outStr.replaceAll(">[^\n]*\n", "").replaceAll("\\s", "");
            final String mapStr = IOUtils.readAll(mutations);
            assertTrue(!genomeTides.equals(stripped));
            final StringBuilder mod = new StringBuilder(stripped);
            final BufferedReader sr = new BufferedReader(new StringReader(mapStr));
            String line;
            int count = 0;
//                      System.err.println("" + genomeTides);
//                      System.err.println("" + xxnums);
//                      System.err.println("" + mod.toString());
            //          boolean insert = true;
            while ((line = sr.readLine()) != null) {
              if (!line.substring(0, 1).equals("#")) {
//                System.err.println(line);
                final String[] parts = line.split("\t");
                assertEquals("sequence1", parts[0]);
                final int pos = Integer.parseInt(parts[1]) - 1;
                assertEquals(".", parts[2]);
                assertEquals("1:3", parts[9]);
                final String refStr = parts[3];
                //final String readStr = parts[4];
                if (refStr.length() < parts[4].length()) {
                  //insert = true;
                  assertEquals(parts[4], mod.toString().substring(pos, pos + parts[4].length()));
                  //              System.err.println("Insert");
//                                System.err.println("pos " + pos + " len " + parts[4].length());
                  mod.delete(pos + 1, pos + parts[4].length());
                } else {
                  //insert = false;
                  // System.err.println("Delete");
                  // System.err.println(readStr);
                  // System.err.println(genomeTides.substring(pos, pos + refStr.length()));
                  assertEquals(refStr, genomeTides.substring(pos, pos + refStr.length()));
                  mod.insert(pos + 1, refStr.substring(1));

                }
                count++;
              }
            }
            assertTrue(stripped.matches("^[ACGT]+$"));
            assertEquals(genomeTides, mod.toString());
            assertTrue(count > 0);
            assertTrue(FileHelper.deleteAll(mutant));
            assertTrue(mutations.delete());
          }
        }
      } finally {
        assertTrue(FileHelper.deleteAll(in));
      }
    } finally {
      assertTrue(FileHelper.deleteAll(temp));
    }
  }

  private boolean snps(final StringBuilder mod, final int pos, final String refBases, final String prediction) {
    final boolean ok = true;
    assertEquals(prediction, mod.substring(pos, pos + prediction.length()));
    mod.delete(pos, pos + prediction.length());
    mod.insert(pos, refBases);
    return ok;
  }

  private boolean deletes(final StringBuilder mod, final int pos, final String refBases, final String prediction) {
    final boolean ok = true;
    assertEquals("i", prediction);
    mod.insert(pos, refBases);
    return ok;
  }

  public void testIndelMutationDiploid() throws IOException {
    final File temp = FileHelper.createTempDirectory();
    try {
      final String genomeTides = "TCTACGCTTAGTACGGCATGCACGATC";
      //final String xxnums      = "012345678901234567890123456";
      final String genome = ">sequence1\n" + genomeTides + "\n";
      final File in = ReaderTestUtils.getDNASubDir(genome, mDir);
      try {
        for (int s = 1; s < 5; s++) {
          try (SequencesReader dsr = SequencesReaderFactory.createDefaultSequencesReader(in)) {
            //System.err.println("===== seed " + s);
            final File mutations = new File(temp, "mutations.mapping");
            final File mutant = new File(temp, "mutant" + s);
            final File twinMutant = new File(temp, "twinmutant");
            final GenomeMutator gm = new GenomeMutator(s, false, false, 1, TEST_PARAMS);
            gm.setRates(0.0, 0.0, 1.0);
            gm.mutateCount(dsr, Sex.EITHER, mutant, twinMutant, new FileOutputStream(mutations), 1);
            //=== mutant
            toFasta(mutant);
            final File output = new File(mutant.getPath() + ".fasta");
            final String outStr = IOUtils.readAll(output);
            final String stripped = outStr.replaceAll(">[^\n]*\n", "").replaceAll("\\s", "");
            final StringBuilder mod = new StringBuilder(stripped);
            //          System.err.println("" + genomeTides);
            //          System.err.println("" + xxnums);
            //          System.err.println("" + mod.toString());
            //=== twin mutant
            toFasta(twinMutant);
            final File output2 = new File(twinMutant.getPath() + ".fasta");
            final String outStr2 = IOUtils.readAll(output2);
            final String stripped2 = outStr2.replaceAll(">[^\n]*\n", "").replaceAll("\\s", "");
            final StringBuilder mod2 = new StringBuilder(stripped2);
            //assertTrue(!(genomeTides.equals(stripped2) && genomeTides.equals(stripped)));

            // test single mutations
            final String mapStr = IOUtils.readAll(mutations);
            final BufferedReader sr = new BufferedReader(new StringReader(mapStr));
            String line;
            int count = 0;
            //          System.err.println("M1:" + genomeTides);
            //          System.err.println("m1:" + xxnums);
            //          System.err.println("M1:" + mod.toString());
            //          System.err.println("M2:" + mod2.toString());
            while ((line = sr.readLine()) != null) {
              if (!line.substring(0, 1).equals("#")) {

                //System.err.println("line:\n" + line);
                final String[] parts = line.split("\t");
                assertEquals(10, parts.length);
                // ref name
                assertEquals("sequence1", parts[0]);
                // pos
                int pos = 0;
                try {
                  pos = Integer.parseInt(parts[1]) - 1;
                } catch (final NumberFormatException e) {
                  fail();
                }
                assertTrue(pos < genomeTides.length() && pos >= 0);
                // type
                assertTrue(parts[9], parts[9].equals("1/1:3") || parts[9].equals("1/0:3"));
                // homozygous do the same twice
                if (parts[3].compareTo("i") != 0) {
                  assertTrue(parts[3].matches("^[ACGT]+$"));
                  assertEquals(parts[3], genomeTides.substring(pos, pos + parts[3].length()));
                  assertTrue(parts[4].matches("^[ACGT]+$") || parts[4].equals("i"));
                  if (!parts[4].equals("i")) {
                    assertTrue(snps(mod, pos, parts[3], parts[4]));
                    assertTrue(snps(mod2, pos, parts[3], parts[4]));
                  } else {
                    assertTrue(deletes(mod, pos, parts[3], parts[4]));
                    assertTrue(deletes(mod2, pos, parts[3], parts[4]));
                  }
                } else {
                  // ref == i
                  //System.errpro
                  assertTrue(parts[4].matches("^[ACGT]+$"));
                  mod.delete(pos, pos + parts[4].length());
                  mod2.delete(pos, pos + parts[4].length());
                }
                count++;
              }
            }
            assertTrue(stripped.matches("^[ACGT]+$"));
            //          System.err.println("X1:" + genomeTides);
            //          System.err.println("x1:" + xxnums);
            //          System.err.println("X1:" + mod.toString());
            //          System.err.println("X2:" + mod2.toString());
            assertEquals(genomeTides, mod.toString());
            assertEquals(genomeTides, mod2.toString());

            assertTrue(stripped2.matches("^[ACGT]+$"));
            // mutator will skip mutations if it has already moved past mutation point
            assertTrue(count > 0);
            assertTrue(FileHelper.deleteAll(mutant));
            assertTrue(FileHelper.deleteAll(twinMutant));
            assertTrue(mutations.delete());
          }
        }
      } finally {
        assertTrue(FileHelper.deleteAll(in));
      }
    } finally {
      assertTrue(FileHelper.deleteAll(temp));
    }
  }

  public void testPriorMutation() throws IOException {
    Diagnostic.setLogStream();
    final File temp = FileHelper.createTempDirectory();
    try {
      final String genome = ""
          + ">a\n"
          + "ACGTACGATCAGCATCTGACATGCTAACGGTCATC" + StringUtils.LS;
      //      + "ACGTACGATCAGCATCTGACATGCTAACGGTCATCGCGGCATTTACGGCACGTCATAGTCGCATGCGGCATTTATAGCGGCGCGTAAATCGT" + "\n"
      //      + ">b\n"
      //      + "CATGATTGCAGCTAGCATCGTGTCACA" + "\n";
      final File in = ReaderTestUtils.getDNASubDir(genome, mDir);
      try {
        try (SequencesReader dsr = SequencesReaderFactory.createDefaultSequencesReader(in)) {
          final GenomeMutator gm = new GenomeMutator(1, true, false, 1, TEST_PARAMS);
          gm.setPriors(mPriors);
          final File mutations = new File(temp, "mutations.mapping");
          final File mutant2 = new File(temp, "mutant2");
          final File twinMutant = new File(temp, "twinmutant");
          final OutputStream fos = new FileOutputStream(mutations);
          try {
            gm.mutatePriors(dsr, Sex.EITHER, mutant2, twinMutant, fos);
          } finally {
            fos.close();
          }
          //String gmstr = gm.toString();
          //System.err.println(gmstr);
          //gen = gm.getMutationGenerator();
        }
      } finally {
        assertTrue(FileHelper.deleteAll(in));
      }
    } finally {
      assertTrue(FileHelper.deleteAll(temp));
    }
  }

  public void testMutationComparator() {
    Diagnostic.setLogStream();
    new MyMutationComparator();
  }

  public void testBurstLimitSuceeds() throws IOException {
    final StringBuilder sb = new StringBuilder(1000);
    sb.append(">name").append(StringUtils.LS);
    for (int i = 0; i < 500; i++) {
      sb.append("N");
    }
    for (int i = 0; i < 125; i++) {
      sb.append("ACGT");
    }
    final File in = ReaderTestUtils.getDNASubDir(sb.toString(), mDir);
    try {
      try (SequencesReader dsr = SequencesReaderFactory.createDefaultSequencesReader(in)) {
        assertEquals(100, countSnpLines(0.1, dsr));
      }
    } finally {
      FileHelper.deleteAll(in);
    }
  }

  public void testBurstLimitFail() throws IOException {
    final StringBuilder sb = new StringBuilder(1000);
    sb.append(">name").append(StringUtils.LS);
    for (int i = 0; i < 800; i++) {
      sb.append("N");
    }
    for (int i = 0; i < 50; i++) {
      sb.append("ACGT");
    }
    final File in = ReaderTestUtils.getDNASubDir(sb.toString(), mDir);
    try {
      try (SequencesReader dsr = SequencesReaderFactory.createDefaultSequencesReader(in)) {
        assertTrue(countSnpLines(0.1, dsr) < 100);
      }
    } finally {
      FileHelper.deleteAll(in);
    }
  }

  private int countSnpLines(double mutationRate, SequencesReader dsr) throws IOException {
    final GenomeMutator gm = new GenomeMutator(1, true, false, 1, TEST_PARAMS);
    final File snps = new File(mDir, "snps");
    final File mutant = new File(mDir, "mutant");
    try {
      try (OutputStream os = new FileOutputStream(snps)) {
        gm.mutateRate(dsr, Sex.EITHER, mutant, null, os, mutationRate);
      }
      return lineCount(snps);
    } finally {
      FileHelper.deleteAll(snps);
      FileHelper.deleteAll(mutant);
    }

  }

  public void testInsertIntoDeleteInvalid() {
    final GenomeMutator gm = new GenomeMutator(1, true, false, 1, TEST_PARAMS);
    final byte[] seq = dnaToSeq("ACCCTTGGGACGTACGTACGTACGT");
    final Mutation delete = new Mutation(1, MutationType.DELETE, false, DifferentMode.HOMOZYGOUS, GenDiffMode.BOTH_SAME, 10, 10, null, null);
    final Mutation insert = new Mutation(3, MutationType.INSERT, false, DifferentMode.HOMOZYGOUS, GenDiffMode.BOTH_SAME, 2, 2, null, null);
    final boolean actual = gm.validPos(new Mutation[] {delete}, 1, insert, seq);
    assertEquals(false, actual);
  }

  public void testInsertIntoMnpInvalid() {
    final GenomeMutator gm = new GenomeMutator(1, true, false, 1, TEST_PARAMS);
    final byte[] seq = dnaToSeq("ACCCTTGGGACGTACGTACGTACGT");
    final Mutation delete = new Mutation(1, MutationType.MNP, true, DifferentMode.DIFFERENT, GenDiffMode.DIFFERENT, 8, 8, null, null);
    final Mutation insert = new Mutation(8, MutationType.INSERT, false, DifferentMode.HOMOZYGOUS, GenDiffMode.BOTH_SAME, 2, 2, null, null);
    final boolean actual = gm.validPos(new Mutation[] {delete}, 1, insert, seq);
    assertEquals(false, actual);

  }

  private void checkArray(String[] exp, String[] actual) {
    assertEquals(exp.length, actual.length);
    for (int i = 0; i < exp.length; i++) {
      assertEquals("exp:" + Arrays.toString(exp) + ",  actual: " + Arrays.toString(actual), exp[i], actual[i]);
    }
  }
  public void testDiffMode() {
    final String[] input = {"Old", "First", "Twin"};
    checkArray(input, GenDiffMode.diffModeArray(input, GenDiffMode.DIFFERENT));
    checkArray(new String[] {"Old", "First"}, GenDiffMode.diffModeArray(input, GenDiffMode.BOTH_SAME));
    checkArray(new String[] {"Old", "First", "Old"}, GenDiffMode.diffModeArray(input, GenDiffMode.FIRST_ONLY));
    checkArray(new String[] {"Old", "Old", "Twin"}, GenDiffMode.diffModeArray(input, GenDiffMode.TWIN_ONLY));
  }

  public void testMutationLineInsertVcf() throws Exception {
    final GenomeMutator gm = new GenomeMutator(1, true, false, 1, new VariantParamsBuilder().create());
    final Mutation m = new Mutation(0, MutationType.INSERT, false, null, null, 1, 1, new byte[] {}, null);
    final VcfRecord rec = gm.mutationLine(null, null, "seq", m, 'N', false, "A", "GA");
    assertEquals("seq\t1\t.\tA\tGA\t.\tPASS\t.\tGT:GQ\t1/1:3", rec.toString());
  }


}
