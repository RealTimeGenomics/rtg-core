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
package com.rtg.reader;

import static com.rtg.util.StringUtils.LS;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.RandomAccessFile;
import java.util.ArrayList;
import java.util.Arrays;

import com.rtg.Slim;
import com.rtg.mode.DNAFastaSymbolTable;
import com.rtg.mode.DnaUtils;
import com.rtg.simulation.genome.RandomDistribution;
import com.rtg.simulation.genome.SequenceGenerator;
import com.rtg.util.PortableRandom;
import com.rtg.util.TestUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.intervals.LongRange;
import com.rtg.util.io.FileUtils;
import com.rtg.util.test.FileHelper;

import junit.framework.Test;
import junit.framework.TestSuite;

/**
 */
public class DefaultSequencesReaderTest extends AbstractSequencesReaderTest {

  @Override
  protected SequencesReader createSequencesReader(final File dir, LongRange region) throws IOException {
    return SequencesReaderFactory.createDefaultSequencesReader(dir, region);
  }

  public static Test suite() {
    return new TestSuite(DefaultSequencesReaderTest.class);
  }

  private File mOutDir = null;

  @Override
  public void setUp() throws Exception {
    super.setUp();
    mOutDir = FileUtils.createTempDir("sequencegenerator", "main");
    Diagnostic.setLogStream();
  }

  @Override
  public void tearDown() throws Exception {
    super.tearDown();
    FileHelper.deleteAll(mOutDir);
    mOutDir = null;
  }

  /**
   * Main to run from tests from command line.
   * @param args ignored.
   */
  public static void main(final String[] args) {
    junit.textui.TestRunner.run(suite());
  }

  private static class TestRaf extends RandomAccessFile {

    private final int mMaxDiff;
    private int mCurrent = 1;
    private int mAdj;

    public TestRaf(final File file, final String mode, final int maxDiff) throws IOException {
      super(file, mode);
      mMaxDiff = maxDiff;
    }

    @Override
    public int read(final byte[] b, final int off, final int len) throws IOException {
      final int len2;
      if (mCurrent < len) {
        len2 = mCurrent;
      } else {
        len2 = len;
      }
      if (mCurrent + mAdj > mMaxDiff) {
        mAdj = -1;
      } else if (mCurrent + mAdj < 1) {
        mAdj = 1;
      }
      mCurrent += mAdj;
      return super.read(b, off, len2);
    }

  }

  public void testLengthsHelper() throws IOException {
    final byte[] buffer = new byte[7];
    final ArrayList<InputStream> al = new ArrayList<>();
    al.add(createStream(">test1\nacgta\n"
                      + ">test2\nagtcatg\n"
                      + ">test3\nacgtttggct\n"
                      + ">test4\natggcttagctacagt\n"
                      + ">test5\nactagattagagtagagatgatgtagatgagtagaaagtt\n"
                      + ">test6\na"));
    //0, 5, 12, 22, 38,
    final FastaSequenceDataSource ds = new FastaSequenceDataSource(al,
            new DNAFastaSymbolTable());
    final SequencesWriter sw = new SequencesWriter(ds, mDir, 20000, PrereadType.UNKNOWN, false);
    sw.processSequences();
    final File file = SdfFileUtils.sequencePointerFile(mDir, 0);
    final int[] exp = {5, 7, 10, 16, 40, 78};
    try (TestRaf raf = new TestRaf(file, "r", 73)) {
      checkRaf(raf, exp, buffer);
    }
    try (RandomAccessFile raf2 = new RandomAccessFile(file, "r")) {
      checkRaf(raf2, exp, buffer);
    }
  }

  private void checkRaf(final RandomAccessFile raf, final int[] expected, final byte[] buffer) throws IOException  {
    final int[] lengths = new int[6];
    DefaultSequencesReader.sequenceLengthsHelper(raf, buffer, lengths, 0, 0, raf.length(), 5);
    assertTrue(Arrays.equals(expected, lengths));
  }

  public void testGenerate() throws Exception {
    final PortableRandom rand = new PortableRandom(1);
    final int[] lens = {2, 5};
    final int[] freq = {1, 1, 1, 1};
    final RandomDistribution rd = new RandomDistribution(freq, rand);
    final SequenceGenerator sdata = new SequenceGenerator(rand, rd, lens, mOutDir);
    final long max = sdata.getSizeLimit();
    assertEquals(1000000000, max);
    sdata.createSequences();
    try (SequencesReader dsr = createSequencesReader(mOutDir)) {
      //System.err.println("SequenceGeneratorTest\n" + dsr.maxLength());
      assertEquals(5, dsr.maxLength());
      assertEquals(2, dsr.minLength());
      assertEquals(2, dsr.numberSequences());

      final String t = "" + dsr.type();
      assertEquals("DNA", t);
      assertEquals("simulatedSequence1", dsr.name(0));
      assertEquals(2, dsr.numberSequences());
      final long len1 = dsr.read(0, new byte[2]);
      assertEquals(2, len1);

      final long len2 = dsr.read(1, new byte[5]);
      assertEquals("simulatedSequence2", dsr.name(1));
      assertEquals(5, len2);
    }

  }

  public void testCopy() throws IOException {
    //set a command line.
    new Slim().intMain(new String[] {"feh", "-f", "super feh"}, TestUtils.getNullOutputStream(), TestUtils.getNullPrintStream());

    final ArrayList<InputStream> al = new ArrayList<>();
    al.add(createStream(">test1\nacgta\n"
                      + ">test2\nagtcatg\n"
                      + ">test3\nacgtttggct\n"
                      + ">test4\natggcttagctacagt\n"
                      + ">test5\nactagattagagtagagatgatgtagatgagtagaaagtt\n"
                      + ">test6\na"));
    //0, 5, 12, 22, 38,
    final FastaSequenceDataSource ds = new FastaSequenceDataSource(al,
            new DNAFastaSymbolTable());
    final SequencesWriter sw = new SequencesWriter(ds, mDir, 20000, PrereadType.UNKNOWN, false);
    sw.setComment("blah rag");
    sw.processSequences();
    try (SequencesReader r = createSequencesReader(mDir)) {
      try (SequencesReader rCopy = r.copy()) {
        assertEquals(r.dataChecksum(), rCopy.dataChecksum());
        assertEquals(r.qualityChecksum(), rCopy.qualityChecksum());
        assertEquals(r.nameChecksum(), rCopy.nameChecksum());
        assertEquals(r.maxLength(), rCopy.maxLength());
        assertEquals(r.minLength(), rCopy.minLength());
        assertEquals(r.numberSequences(), rCopy.numberSequences());
        assertEquals(r.sdfVersion(), rCopy.sdfVersion());
        if (r instanceof AnnotatedSequencesReader) {
          final AnnotatedSequencesReader rCast = (AnnotatedSequencesReader) r;
          final AnnotatedSequencesReader rCopyCast = (AnnotatedSequencesReader) rCopy;
          assertEquals(rCast.comment(), rCopyCast.comment());
          assertEquals("feh -f \"super feh\"", rCopyCast.commandLine());
        }
      }
    }

  }

  static final String FASTA = ">r0" + LS + "ACGTACG" + LS
   + ">r1" + LS + "GCGTA" + LS
   + ">r2" + LS + "cCGTAC" + LS
   + ">r3" + LS + "TCGTACGTAC" + LS
   + ">r4" + LS + "GGGTACGTACGT" + LS;
  public void testRegion() throws IOException {
    ReaderTestUtils.getReaderDNA(FASTA, mOutDir, new SdfId(0L)).close();
    try (DefaultSequencesReader reader = new DefaultSequencesReader(mOutDir, new LongRange(1, 3))) {
      assertEquals(2, reader.numberSequences());
      byte[] data = reader.read(0);
      assertEquals("GCGTA", DnaUtils.bytesToSequenceIncCG(data));
      data = reader.read(1);
      assertEquals("CCGTAC", DnaUtils.bytesToSequenceIncCG(data));

      final int[] lengths = reader.sequenceLengths(0, 2);
      assertEquals(2, lengths.length);
      assertEquals(5, lengths[0]);
      assertEquals(6, lengths[1]);
    }


  }


  public void testReadMe() throws IOException {
    ReaderTestUtils.getReaderDNA(FASTA, mOutDir, new SdfId(0L)).close();
    try (DefaultSequencesReader reader = new DefaultSequencesReader(mOutDir, new LongRange(1, 3))) {
      assertNull(reader.getReadMe());
      FileUtils.stringToFile("This is some example content\nBlah", new File(mOutDir, "readme.txt"));
      assertEquals("This is some example content\nBlah", reader.getReadMe());
    }
  }

  public void testSeekNegative() throws IOException {
    ReaderTestUtils.getReaderDNA(FASTA, mOutDir, new SdfId(0L)).close();
    try (final DefaultSequencesReader reader = new DefaultSequencesReader(mOutDir, new LongRange(1, 3))) {
      try {
        reader.iterator().seek(Integer.MIN_VALUE);
        fail();
      } catch (final IllegalArgumentException e) {
        // ok
      }
    }
  }

}

