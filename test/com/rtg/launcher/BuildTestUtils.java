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
package com.rtg.launcher;


import java.io.ByteArrayInputStream;
import java.io.File;
import java.io.FilenameFilter;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;

import com.rtg.mode.DNAFastaSymbolTable;
import com.rtg.reader.FastaSequenceDataSource;
import com.rtg.reader.PrereadType;
import com.rtg.reader.SequencesReader;
import com.rtg.reader.SequencesReaderFactory;
import com.rtg.reader.SequencesWriter;
import com.rtg.util.StringUtils;
import com.rtg.util.test.FileHelper;

import junit.framework.TestCase;

/**
 * Utilities for helping the tests of the various build tasks.
 */
public final class BuildTestUtils {

  private BuildTestUtils() { }

  /**
   * @param inputSequence the sequence
   * @param dir the dir
   * @throws IOException if badness
   */
  public static void writeDNA(final String inputSequence, final File dir) throws IOException {
    final ArrayList<InputStream> inputStreams = new ArrayList<>();
    inputStreams.add(new ByteArrayInputStream(inputSequence.getBytes()));
    final FastaSequenceDataSource ds = new FastaSequenceDataSource(inputStreams,
        new DNAFastaSymbolTable());
    final SequencesWriter sequenceWriter = new SequencesWriter(ds, dir, 20, PrereadType.UNKNOWN, true);
    sequenceWriter.processSequences();
  }

  /**
   * @param inputSequence the sequence
   * @param dir the dir
   * @return the reader
   * @throws IOException if badness
   */
  public static SequencesReader getReaderDNA(final String inputSequence, final File dir) throws IOException {
    writeDNA(inputSequence, dir);
    return SequencesReaderFactory.createDefaultSequencesReader(dir);
  }

  /**
   * Make sure there are no logs in the working directory.
   */
  public static void checkNoLogs() {
    final File wd = new File(System.getProperty("user.dir"));
    final File[] logs = wd.listFiles(new FilenameFilter() {
      @Override
      public boolean accept(final File dir, final String name) {
        return name.startsWith("log") && name.endsWith(".log");
      }
    });
    if (logs != null && logs.length > 0) {
      System.err.println("Found " + logs.length + " logs:");
      for (final File log : logs) {
        System.err.println(log);
      }
      TestCase.fail();
    }
  }

  /**
   * Take an array of strings which may contain nulls and return one with
   * the nulls elided.
   * Useful when generating command line arguments.
   * @param args0 to be processed.
   * @return args0 without the nulls.
   */
  public static String[] elideNulls(final String[] args0) {
    int count = 0;
    for (final String str : args0) {
      if (str != null) {
        count++;
      }
    }
    final String[] args = new String[count];
    int i = 0;
    for (final String str : args0) {
      if (str != null) {
        args[i] = str;
        i++;
      }
    }
    return args;
  }

  //methods stolen from BuildSearchUtils in internal
  /** Subject sequence used for the calibration runs.  */
  public static final String SEQ_DNA_A = ""
    + ">x" + StringUtils.LS
    + "actg" + StringUtils.LS;
  /**
   * Pre-read a simple sequence and return the (temporary) directory where it has been written.
   * @param parent the parent dir
   * @return the directory where the sequence is stored.
   * @throws IOException if an I/O error occurs.
   */
  public static File prereadDNA(File parent) throws IOException {
    final File dir = FileHelper.createTempDirectory(parent);
    try {
      prereadDNA(SEQ_DNA_A, dir);
    } catch (IOException e) {
      FileHelper.deleteAll(dir);
      throw e;
    }
    return dir;
  }

  /**
   * Pre-read a sequence and return the (temporary) directory where it has been written.
   * @param parent the parent dir
   * @param inputDnaSequence data in fasta format.
   * @return the directory where the sequence is stored.
   * @throws IOException if an I/O error occurs.
   */
  public static File prereadDNA(File parent, final String inputDnaSequence) throws IOException {
    final File dir = FileHelper.createTempDirectory(parent);
    try {
      prereadDNA(inputDnaSequence, dir);
    } catch (IOException e) {
      FileHelper.deleteAll(dir);
      throw e;
    }
    return dir;
  }


  /**
   * Pre-read a sequence and return the(temporary) directory where it has been written.
   * @param inputDnaSequence data in fasta format.
   * @param dir where to write the pre-read data.
   * @throws IOException if an I/O error occurs.
   */
  public static void prereadDNA(final String inputDnaSequence, final File dir) throws IOException {
    final ArrayList<InputStream> inputStreams = new ArrayList<>();
    inputStreams.add(new ByteArrayInputStream(inputDnaSequence.getBytes()));
    final FastaSequenceDataSource ds =
      new FastaSequenceDataSource(
          inputStreams,
          new DNAFastaSymbolTable()
      );
    final SequencesWriter sequenceWriter = new SequencesWriter(ds, dir, 20, PrereadType.UNKNOWN, true);
    sequenceWriter.processSequences();
  }

}
