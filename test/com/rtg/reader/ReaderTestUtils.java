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

import java.io.ByteArrayInputStream;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.UUID;

import com.rtg.mode.DNAFastaSymbolTable;
import com.rtg.mode.ProteinFastaSymbolTable;
import com.rtg.reader.FastqSequenceDataSource.FastQScoreType;
import com.rtg.util.Constants;
import com.rtg.util.intervals.LongRange;
import com.rtg.util.StringUtils;
import com.rtg.util.io.FileUtils;


/**
 * Reader-based utilities for when testing.
 */
public final class ReaderTestUtils {

  /**
   * The dummy SDF ID that will be put into SDF files created by this,
   * to be used in the construction of files so that they will match.
   */
  public static final SdfId DUMMY_TEST_ID = new SdfId(new UUID(0, 42));

  /** Prevent construction. */
  private ReaderTestUtils() {
  }

  /**
   * Construct a pre-reader directory containing the specified fasta sequence.
   * @param inputDnaSequence the sequence.
   * @param dir directory where pre-reader placed.
   * @param sdfId Id to use for SDF, <code>null</code> for default
   * @return SequencesReader using the directory.
   * @throws IOException whenever.
   */
  public static SequencesReader getReaderDNA(final String inputDnaSequence, final File dir, SdfId sdfId) throws IOException {
    return getReaderDNA(inputDnaSequence, dir, sdfId, false);
  }

  /**
   * Construct a pre-reader directory containing the specified fasta sequence.
   * @param inputDnaSequenceLeft the left sequence.
   * @param inputDnaSequenceRight the right sequence.
   * @param dir SDF dir
   * @param sdfId Id to use for SDF, <code>null</code> for default
   * @return SDF-ID actually used for SDF
   * @throws IOException if an IO Error occured
   */
  public static SdfId createPairedReaderDNA(final String inputDnaSequenceLeft, final String inputDnaSequenceRight, final File dir, SdfId sdfId) throws IOException {
    final SdfId id;
    try (SequencesReader sr1 = getReaderDNA(inputDnaSequenceLeft, new File(dir, "left"), sdfId)) {
      id = sr1.getSdfId();
      final SequencesReader sr2 = getReaderDNA(inputDnaSequenceRight, new File(dir, "right"), id);
      sr2.close();
    }
    return id;
  }

  /**
   * Construct a pre-reader directory containing the specified fasta sequence.
   * @param inputDnaSequence the sequence.
   * @param dir directory where pre-reader placed.
   * @param sdfId Id to use for SDF, <code>null</code> for default
   * @param memSeqReader use memory sequences Reader
   * @return SequencesReader using the directory.
   * @throws IOException whenever.
   */
  public static SequencesReader getReaderDNA(final String inputDnaSequence, final File dir, SdfId sdfId, final boolean memSeqReader) throws IOException {
    final ArrayList<InputStream> inputStreams = new ArrayList<>();
    inputStreams.add(new ByteArrayInputStream(inputDnaSequence.getBytes()));
    final FastaSequenceDataSource ds = new FastaSequenceDataSource(inputStreams, new DNAFastaSymbolTable());
    final SequencesWriter sequenceWriter = new SequencesWriter(ds, dir, Constants.MAX_FILE_SIZE, PrereadType.UNKNOWN, true);
    if (sdfId != null) {
      sequenceWriter.setSdfId(sdfId);
    }
    sequenceWriter.processSequences();
    return memSeqReader ? SequencesReaderFactory.createMemorySequencesReader(dir, true, LongRange.NONE) : SequencesReaderFactory.createDefaultSequencesReader(dir);
  }

  /**
   * Construct a pre-reader directory containing the specified fasta sequence.
   * @param inputDnaSequence the sequence.
   * @param dir directory where pre-reader placed.
   * @param sdfId Id to use for SDF, <code>null</code> for default
   * @param maxFileSize maximum file size for each chunk
   * @return SequencesReader using the directory.
   * @throws IOException whenever.
   */
  public static SequencesReader getReaderDNA(final String inputDnaSequence, final File dir, SdfId sdfId, final long maxFileSize) throws IOException {
    final ArrayList<InputStream> inputStreams = new ArrayList<>();
    inputStreams.add(new ByteArrayInputStream(inputDnaSequence.getBytes()));
    final FastaSequenceDataSource ds = new FastaSequenceDataSource(inputStreams, new DNAFastaSymbolTable());
    final SequencesWriter sequenceWriter = new SequencesWriter(ds, dir, maxFileSize, PrereadType.UNKNOWN, true);
    if (sdfId != null) {
      sequenceWriter.setSdfId(sdfId);
    }
    sequenceWriter.processSequences();
    return SequencesReaderFactory.createDefaultSequencesReader(dir);
  }

    /**
   * Construct a pre-reader directory containing the specified fasta sequence.
   * @param inputDnaSequence the sequence.
   * @param dir directory where pre-reader placed.
   * @param sdfId Id to use for SDF, <code>null</code> for default
   * @param maxFileSize maximum file size for each chunk
   * @return SequencesReader using the directory.
   * @throws IOException whenever.
   */
  public static SequencesReader getReaderDNAFastq(final String inputDnaSequence, final File dir, SdfId sdfId, final long maxFileSize) throws IOException {
    final ArrayList<InputStream> inputStreams = new ArrayList<>();
    inputStreams.add(new ByteArrayInputStream(inputDnaSequence.getBytes()));
    final FastaSequenceDataSource ds = new FastqSequenceDataSource(inputStreams, FastQScoreType.PHRED);
    final SequencesWriter sequenceWriter = new SequencesWriter(ds, dir, maxFileSize, PrereadType.UNKNOWN, true);
    if (sdfId != null) {
      sequenceWriter.setSdfId(sdfId);
    }
    sequenceWriter.processSequences();
    return SequencesReaderFactory.createDefaultSequencesReader(dir);
  }

  /**
   * Construct a pre-reader directory containing the specified fasta sequence.
   * @param inputDnaSequence the sequence.
   * @param dir directory where pre-reader placed.
   * @param memSeqReader use memory sequences Reader
   * @param guid the guid for this input
   * @return SequencesReader using the directory.
   * @throws IOException whenever.
   */
  public static SequencesReader getReaderDNA(final String inputDnaSequence, final File dir, final boolean memSeqReader, final SdfId guid) throws IOException {
    final ArrayList<InputStream> inputStreams = new ArrayList<>();
    inputStreams.add(new ByteArrayInputStream(inputDnaSequence.getBytes()));
    final FastaSequenceDataSource ds = new FastaSequenceDataSource(inputStreams, new DNAFastaSymbolTable());
    final SequencesWriter sequenceWriter = new SequencesWriter(ds, dir, Constants.MAX_FILE_SIZE, PrereadType.UNKNOWN, true);
    sequenceWriter.setSdfId(guid);
    sequenceWriter.processSequences();
    return memSeqReader ? SequencesReaderFactory.createMemorySequencesReader(dir, true, LongRange.NONE) : SequencesReaderFactory.createDefaultSequencesReader(dir);
  }

  /**
   * @param sequence fasta sequence
   * @return sequences reader for given sequence
   * @throws IOException if data is invalid
   */
  public static SequencesReader getReaderDnaMemory(String sequence) throws IOException {
    final SequencesWriter sw = new SequencesWriter(new FastaSequenceDataSource(Arrays.asList((InputStream) new ByteArrayInputStream(sequence.getBytes())), new DNAFastaSymbolTable()), null, PrereadType.UNKNOWN, true);
    sw.setSdfId(DUMMY_TEST_ID);
    return sw.processSequencesInMemory(null, true, new SimplePrereadNames(), new SimplePrereadNames(), LongRange.NONE);
  }

  /**
   * Produce a fasta string with the given sequences
   * @param sequences a list of sequences
   * @return a fasta string with all the sequences
   */
  public static String fasta(String... sequences) {
    final StringBuilder sb = new StringBuilder();
    for (int i = 0; i < sequences.length; i++) {
      sb.append(">").append(i).append(StringUtils.LS).append(sequences[i]).append(StringUtils.LS);
    }
    return sb.toString();
  }

  /**
   * Construct a pre-reader directory containing the specified fastq sequence.
   * @param inputDnaSequence the sequence.
   * @param dir directory where pre-reader placed.
   * @param isSolexa is Solexa
   * @return SequencesReader using the directory.
   * @throws IOException whenever.
   */
  public static SequencesReader getReaderDNAFastq(final String inputDnaSequence, final File dir, final boolean isSolexa) throws IOException {
    return getReaderDNAFastq(inputDnaSequence, dir, Constants.MAX_FILE_SIZE, isSolexa);
  }

  /**
   * Construct a pre-reader directory containing the specified fastq sequence.
   * @param inputDnaSequence the sequence.
   * @param dir directory where pre-reader placed.
   * @param maxFileSize maximum size for output files
   * @param isSolexa is Solexa
   * @return SequencesReader using the directory.
   * @throws IOException whenever.
   */
  public static SequencesReader getReaderDNAFastq(final String inputDnaSequence, final File dir, final long maxFileSize, final boolean isSolexa) throws IOException {
    final ArrayList<InputStream> inputStreams = new ArrayList<>();
    inputStreams.add(new ByteArrayInputStream(inputDnaSequence.getBytes()));
    final FastqSequenceDataSource ds = new FastqSequenceDataSource(inputStreams, isSolexa ? FastQScoreType.SOLEXA : FastQScoreType.PHRED);
    final SequencesWriter sequenceWriter = new SequencesWriter(ds, dir, maxFileSize, PrereadType.UNKNOWN, true);
    sequenceWriter.setSdfId(DUMMY_TEST_ID);
    sequenceWriter.processSequences();
    return SequencesReaderFactory.createDefaultSequencesReader(dir);
  }


  /**
   * Construct an SDF directory containing the specified fastq sequence as CG data.
   * @param inputDnaSequence the sequence.
   * @param dir directory where SDF placed.
   * @param prereadArm some special sauce.
   * @return SequencesReader using the directory.
   * @throws IOException whenever.
   */
  public static SequencesReader getReaderDNAFastqCG(final String inputDnaSequence, final File dir, final PrereadArm prereadArm) throws IOException {
    final ArrayList<InputStream> inputStreams = new ArrayList<>();
    inputStreams.add(new ByteArrayInputStream(inputDnaSequence.getBytes()));
    final FastqSequenceDataSource ds = new FastqSequenceDataSource(inputStreams, FastQScoreType.PHRED);
    final SequencesWriter sequenceWriter = new SequencesWriter(ds, dir, Constants.MAX_FILE_SIZE, PrereadType.CG, true);
    sequenceWriter.setSdfId(DUMMY_TEST_ID);
    sequenceWriter.setPrereadArm(prereadArm);
    sequenceWriter.processSequences();
    return SequencesReaderFactory.createDefaultSequencesReader(dir);
  }

  /**
   * Construct an SDF directory containing the specified sequence.
   * @param inputProteinSequence the sequence.
   * @param dir directory where SDF placed.
   * @return SequencesReader using the directory.
   * @throws IOException whenever.
   */
  public static SequencesReader getReaderProtein(final String inputProteinSequence, final File dir) throws IOException {

    final ArrayList<InputStream> inputStreams = new ArrayList<>();
    inputStreams.add(new ByteArrayInputStream(inputProteinSequence.getBytes()));
    final FastaSequenceDataSource ds = new FastaSequenceDataSource(inputStreams, new ProteinFastaSymbolTable());
    final SequencesWriter sequenceWriter = new SequencesWriter(ds, dir, 100000, PrereadType.UNKNOWN, true);
    sequenceWriter.setSdfId(DUMMY_TEST_ID);
    sequenceWriter.processSequences();

    return SequencesReaderFactory.createDefaultSequencesReader(dir);
  }

  /**
   * Construct a temporary SDF directory containing the specified sequence.
   * Note: it is the callers responsibility to delete the directory ( FileUtils.deleteFiles is a possibility).
   * @param inputDnaSequence the sequence.
   * @return SequencesReader using the directory.
   * @throws IOException whenever.
   */
  public static File getDNADir(final String inputDnaSequence) throws IOException {
    final File dir = FileUtils.createTempDir("testfasta", "dna");
    return getDNADir(inputDnaSequence, dir);
  }
  /**
   * Construct a temporary SDF directory containing the specified sequence.
   * Note: it is the callers responsibility to delete the directory ( FileUtils.deleteAll is a possibility).
   * @param inputDnaSequence the sequence.
   * @param dir the directory to put the SDF in.
   * @return SequencesReader using the directory.
   * @throws IOException whenever.
   */
  public static File getDNADir(final String inputDnaSequence, File dir) throws IOException {
    return getDNADir(inputDnaSequence, dir, true, true, false);
  }

  /**
   * Construct a temporary SDF directory containing the specified sequence.
   * Note: it is the callers responsibility to delete the directory ( FileUtils.deleteAll is a possibility).
   * @param inputDnaSequence the sequence.
   * @param dir the directory to put the SDF in.
   * @param includeQual include quality values
   * @param includeName include name values
   * @param includeRef include reference file
   * @return SequencesReader using the directory.
   * @throws IOException whenever.
   */
  public static File getDNADir(final String inputDnaSequence, File dir, boolean includeQual, boolean includeName, boolean includeRef) throws IOException {
    final ArrayList<InputStream> inputStreams = new ArrayList<>();
    inputStreams.add(new ByteArrayInputStream(inputDnaSequence.getBytes()));
    final FastaSequenceDataSource ds = new FastaSequenceDataSource(inputStreams, new DNAFastaSymbolTable());
    final SequencesWriter sequenceWriter = new SequencesWriter(ds, dir, Constants.MAX_FILE_SIZE, PrereadType.UNKNOWN, true);
    sequenceWriter.setSdfId(DUMMY_TEST_ID);
    sequenceWriter.processSequences(includeQual, includeName);
    if (includeRef) {
      final File referenceFile = new File(dir, "reference.txt");
      try (PrintStream fw = new PrintStream(referenceFile)) {
        fw.println("version 0");
        fw.println("either def diploid linear");
      }
    }
    return dir;
  }

  /**
   * Construct a temporary sdf directory containing the specified sequence.
   * Note: it is the callers repsonsibility to delete the directory ( FileUtils.deleteAll is a possibility).
   * @param inputDnaSequence the sequence.
   * @param isSolexa if Solexa input
   * @return SequencesReader using the directory.
   * @throws IOException whenever.
   */
  public static File getDNAFastqDir(final String inputDnaSequence, final boolean isSolexa) throws IOException {
    final File dir = FileUtils.createTempDir("testfastq", "dna");
    return getDNAFastqDir(inputDnaSequence, dir, isSolexa);
  }

  /**
   * Construct a temporary sdf directory containing the specified sequence.
   * Note: it is the callers repsonsibility to delete the directory ( FileUtils.deleteAll is a possibility).
   * @param inputDnaSequence the sequence.
   * @param dir directory to put sdf in
   * @param isSolexa if Solexa input
   * @return SequencesReader using the directory.
   * @throws IOException whenever.
   */
  public static File getDNAFastqDir(final String inputDnaSequence, final File dir, final boolean isSolexa) throws IOException {
    final ArrayList<InputStream> inputStreams = new ArrayList<>();
    inputStreams.add(new ByteArrayInputStream(inputDnaSequence.getBytes()));
    final FastqSequenceDataSource ds = new FastqSequenceDataSource(inputStreams, isSolexa ? FastQScoreType.SOLEXA : FastQScoreType.PHRED);
    final SequencesWriter sequenceWriter = new SequencesWriter(ds, dir, Constants.MAX_FILE_SIZE, PrereadType.UNKNOWN, true);
    sequenceWriter.setSdfId(DUMMY_TEST_ID);
    sequenceWriter.processSequences();
    return dir;
  }
}
