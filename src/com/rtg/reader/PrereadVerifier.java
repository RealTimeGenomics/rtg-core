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

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.util.Arrays;

import com.rtg.launcher.CommonFlags;
import com.rtg.util.cli.CFlags;
import com.rtg.util.cli.Validator;
import com.rtg.util.diagnostic.CliDiagnosticListener;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.ErrorType;
import com.rtg.util.diagnostic.SlimException;
import com.rtg.util.intervals.LongRange;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.LogFile;
import com.rtg.util.io.LogStream;

/**
 */
public final class PrereadVerifier {

  private PrereadVerifier() {
  }

  private static final Validator VALIDATOR = new Validator() {
    @Override
    public boolean isValid(final CFlags flags) {
      final File input = (File) flags.getAnonymousValue(0);
      return CommonFlags.validateSDF(input);
    }
  };

  private static final String APPLICATION_NAME = "SDFVerify";

  /**
   * @param args command line arguments
   */
  public static void main(final String[] args) {
    System.exit(mainInit(args, System.out, System.err));
  }

  /**
   * Deal with initial log file.
   *
   * @param args command line arguments
   * @param out stdout where Duster output is written.
   * @param err where error and warning messages are written for user consumption.
   * @return exit code - 0 if all ok - 1 if command line arguments failed.
   */
  static int mainInit(final String[] args, final PrintStream out, final PrintStream err) {
    return mainExec(args, out, err, PrereadVerifier.initialLog(APPLICATION_NAME));
  }

  static int mainExec(final String[] args, final PrintStream out, final PrintStream err, final LogStream initLog) {
    CliDiagnosticListener listener = null;
    try {
      listener = PrereadVerifier.initializeLogs(args, err, initLog);
      final CFlags flags = getCFlags(APPLICATION_NAME, out, err);
      if (flags.setFlags(args)) {
        if (mainParsed(flags, out)) {
          Diagnostic.deleteLog(); // was successful execution
          return 0;
        } else {
          return 1;
        }
      } else {
        Diagnostic.deleteLog(); // was command line error
        return 1;
      }
    } catch (final SlimException e) {
      throw e;
    } catch (final Throwable t) { // catch everything except SlimException
      throw new SlimException(t);
    } finally { // make sure listener removed so later unit tests not compromised
      Diagnostic.removeListener(listener);
    }
  }

  static boolean mainParsed(final CFlags flags, final PrintStream out) throws IOException {
    final File input = (File) flags.getAnonymousValue(0);
    //final boolean fast = flags.isSet(FAST);
    final File inleft = new File(input, "left");
    final File inright = new File(input, "right");
    if (inleft.exists() && inright.exists()) {
      if (verifyDir(inleft) && verifyDir(inright)) {
        out.println("\nPaired-end SDF verified okay.");
        return true;
      } else {
        return false;
      }
    } else if (verifyDir(input)) {
      out.println("\nSingle-end SDF verified okay.");
      return true;
    } else {
      return false;
    }
  }

  /**
   * Method to verify <code>Preread</code> directory
   * @param f file
   * @return true if all okay
   * @throws IOException if an IO error occurs
   */
  public static boolean verifyDir(final File f) throws IOException {
    try (SequencesReader reader = SequencesReaderFactory.createMemorySequencesReader(f, true, true, LongRange.NONE)) {
      //reader = SequencesReaderFactory.createDefaultSequencesReader(f);
      //reader.globalIntegrity();
      return verify(reader, f);
    } catch (final IOException e) {
      Diagnostic.error(ErrorType.SDF_VERIFICATION_FAILED);
      return false;
    } catch (final RuntimeException ex) {
      //failed to create reader
      Diagnostic.error(ErrorType.SDF_VERIFICATION_FAILED);
      return false;
    }
  }

  static boolean verify(final SequencesReader reader, final File input) {
    if (reader.sdfVersion() < IndexFile.SINGLE_CHECKSUM_VERSION) {
      Diagnostic.error(ErrorType.SDF_VERSION_INVALID, input.getAbsolutePath());
      return false;
    }
    try {
      if (reader.sdfVersion() < IndexFile.SEPARATE_CHECKSUM_VERSION) {
        //System.out.println("Verifying version 2-6 checksum");
        final long checksum = reader.dataChecksum();
        final long newChecksum = calcChecksumOld(reader);
        if (newChecksum == 0) {
          return false;
        }
        if (newChecksum != checksum) {
          Diagnostic.error(ErrorType.SDF_VERIFICATION_FAILED);
          return false;
        }
      } else {
        //System.out.println("Verifying version 7+ data checksum");
        long checksum = reader.dataChecksum();
        long newChecksum = calcDataChecksum(reader);
        if (newChecksum == 0) {
          return false;
        }
        if (newChecksum != checksum) {
          Diagnostic.error(ErrorType.SDF_VERIFICATION_FAILED);
          return false;
        }

        if (reader.hasQualityData()) {
          //System.out.println("Verifying version 7+ quality checksum");
          checksum = reader.qualityChecksum();
          newChecksum = calcQualityChecksum(reader);
          if (newChecksum == 0) {
            return false;
          }
          if (newChecksum != checksum) {
            Diagnostic.error(ErrorType.SDF_VERIFICATION_FAILED);
            return false;
          }
        }

        //System.out.println("Verifying version 7+ names checksum");
        checksum = reader.nameChecksum();
        newChecksum = calcNameChecksum(reader);
        if (newChecksum == 0) {
          return false;
        }
        if (newChecksum != checksum) {
          Diagnostic.error(ErrorType.SDF_VERIFICATION_FAILED);
          return false;
        }
        //check suffix checksum
        checksum = reader.suffixChecksum();
        newChecksum = calcNameSuffixChecksum(reader);
        if (newChecksum == 0) {
          return false;
        }
        if (newChecksum != checksum) {
          Diagnostic.error(ErrorType.SDF_VERIFICATION_FAILED);
          return false;
        }
      }
    } catch (final Throwable e) {
      Diagnostic.error(ErrorType.SDF_VERIFICATION_FAILED);
      return false;
    }
    return true;
  }

  private static long calcChecksumOld(final SequencesReader reader) throws IllegalStateException, IOException {
    final PrereadHashFunction prf = new PrereadHashFunction();
    long totalDone = 0;
    final long totalNumberOfSeq = reader.numberSequences();
    final long totalTides = reader.totalLength();
    for (long numberOfSeq = 0; numberOfSeq < totalNumberOfSeq; numberOfSeq++) {
      if (numberOfSeq >= totalNumberOfSeq) {
        return 0;
      }
      final int currentSeqLen = reader.length(numberOfSeq);
      // Safety for OOM condition in case of corrupt length
      if (currentSeqLen < 0 || currentSeqLen > totalTides) {
        return 0;
      }

      final byte[] data = new byte[currentSeqLen];
      final int size = reader.read(numberOfSeq, data);
      if (size != currentSeqLen) {
        return 0;
      }
      totalDone += size;
      if (totalDone > totalTides) {
        return 0;
      }
      for (final byte b : data) {
        prf.irvineHash(b);
      }
      if (reader.hasQualityData()) {
        reader.readQuality(numberOfSeq, data);
        if (size != currentSeqLen) {
          return 0;
        }
        for (final byte b : data) {
          prf.irvineHash(b);
        }
      }
      prf.irvineHash(reader.name(numberOfSeq));
    }
    try (FileInputStream seqIndexIn = new FileInputStream(SdfFileUtils.sequenceIndexFile(reader.path()))) {
      prf.irvineHash(seqIndexIn);
    }
    try (FileInputStream labIndexIn = new FileInputStream(SdfFileUtils.labelIndexFile(reader.path()))) {
      prf.irvineHash(labIndexIn);
    }
    return prf.getHash();
  }


  private static long calcDataChecksum(final SequencesReader reader) throws IOException {
    final PrereadHashFunction dataf = new PrereadHashFunction();
    final long totalNumberOfSeq = reader.numberSequences();
    final long totalTides = reader.totalLength();
    long totalDone = 0;
    for (long numberOfSeq = 0; numberOfSeq < totalNumberOfSeq; numberOfSeq++) {
      if (numberOfSeq >= totalNumberOfSeq) {
        return 0;
      }
      final int currentSeqLen = reader.length(numberOfSeq);
      // Safety for OOM condition in case of corrupt length
      if (currentSeqLen < 0 || currentSeqLen > totalTides) {
        return 0;
      }

      final byte[] data = new byte[currentSeqLen];
      final int size = reader.read(numberOfSeq, data);
      if (size != currentSeqLen) {
        return 0;
      }
      totalDone += size;
      if (totalDone > totalTides) {
        return 0;
      }
      for (final byte b : data) {
        dataf.irvineHash(b);
      }
      dataf.irvineHash((long) currentSeqLen);
    }
    return dataf.getHash();
  }

  private static long calcQualityChecksum(final SequencesReader reader) throws IOException {
    final PrereadHashFunction qualf = new PrereadHashFunction();
    final long totalNumberOfSeq = reader.numberSequences();
    final long totalTides = reader.totalLength();
    long totalDone = 0;
    for (long numberOfSeq = 0; numberOfSeq < totalNumberOfSeq; numberOfSeq++) {
      if (numberOfSeq >= totalNumberOfSeq) {
        return 0;
      }
      final int currentSeqLen = reader.length(numberOfSeq);
      // Safety for OOM condition in case of corrupt length
      if (currentSeqLen < 0 || currentSeqLen > totalTides) {
        return 0;
      }

      final byte[] data = new byte[currentSeqLen];
      final int size = reader.readQuality(numberOfSeq, data);
      if (size != currentSeqLen) {
        return 0;
      }
      totalDone += size;
      if (totalDone > totalTides) {
        return 0;
      }
      for (final byte b : data) {
        qualf.irvineHash(b);
      }
      qualf.irvineHash((long) currentSeqLen);
    }
    return qualf.getHash();
  }

  private static long calcNameChecksum(final SequencesReader reader) throws IOException {
    final PrereadHashFunction namef = new PrereadHashFunction();
    final long totalNumberOfSeq = reader.numberSequences();
    long numberOfSeq = -1;
    final SequencesIterator it = reader.iterator();
    it.seek(0);
    do {
      numberOfSeq++;
      if (numberOfSeq >= totalNumberOfSeq) {
        return 0;
      }
      namef.irvineHash(it.currentName());
      namef.irvineHash(it.currentName().length());
    } while (it.nextSequence());
    return namef.getHash();
  }

  private static long calcNameSuffixChecksum(final SequencesReader reader) throws IOException {
    final PrereadHashFunction namef = new PrereadHashFunction();
    final long totalNumberOfSeq = reader.numberSequences();
    long numberOfSeq = -1;
    final SequencesIterator it = reader.iterator();
    it.seek(0);
    do {
      numberOfSeq++;
      if (numberOfSeq >= totalNumberOfSeq) {
        return 0;
      }
      namef.irvineHash(it.currentNameSuffix());
      namef.irvineHash(it.currentNameSuffix().length());
    } while (it.nextSequence());
    return namef.getHash();
  }


  private static CFlags getCFlags(final String name, final Appendable out, final Appendable err) {
    final CFlags flags = new CFlags(name, out, err);
    flags.registerRequired(File.class, "SDF", "the SDF to be verified");
    flags.setValidator(VALIDATOR);
    return flags;
  }

  /**
   * Do common tasks for initializing logs.
   * @param args command line arguments
   * @param err errors and warnings for user, also final last resort place to put logs if cannot create one in working directory.
   * @param initLog where to send the inital logging.
   * @return the listener that writes the errors and warnings for the user.
   */
  public static CliDiagnosticListener initializeLogs(final String[] args, final PrintStream err, final LogStream initLog) {
    if (initLog != null) {
      Diagnostic.setLogStream(initLog);
      Diagnostic.logEnvironment();
    } else {
      Diagnostic.setLogStream();
    }
    Diagnostic.userLog("Command line arguments:" + Arrays.toString(args));
    final CliDiagnosticListener listener = new CliDiagnosticListener(err);
    Diagnostic.addListener(listener);
    return listener;
  }

  /**
   * Create the place to write the initial log before any command line arguments are parsed (and we know where to put the real logs).
   * @param name name of application requesting logging
   * @return the <code>LogStream</code> to use for the initial logging.
   */
  public static LogStream initialLog(final String name) {
    try {
      return new LogFile(File.createTempFile(name.replace(' ', '-'), "." + FileUtils.LOG_SUFFIX, new File(System.getProperty("user.dir"))));
    } catch (final IOException e) {
      return null;
    }
  }
}
